#include "prims/KDPrims.h"

namespace scds {

constexpr int num_bins = 1;
constexpr int max_depth_allowed   = 1;
constexpr int split_node_prim_num = 8;
constexpr float traverse_cost = 1;

template <typename Ty>
void recursive_kd_tree_SAH(KDPrimNode<Ty>* const cur_node, std::vector<AABB<Ty>>& aabb_infos, int depth = 0) {
    // tree depth constrains and minimum primitive num constrains
    if (cur_node->num_elems < split_node_prim_num || depth > max_depth_allowed) return;
    AABB<Ty> fwd_bound, bwd_bound;
    int child_prim_cnt = 0, seg_bin_idx = 0;
    Ty min_cost  = 5e9, node_prim_cnt = Ty(cur_node->num_elems), 
       split_pos = 0, node_inv_area = 1. / cur_node->bound.area();
    SplitAxis max_axis = NONE;

    {   // local scope, destroy bins and idx_bins before the recursive calls
        std::vector<Ty> bins;
        std::vector<AxisBins<Ty>> idx_bins;
        // Step 1: decide the axis that expands the maximum extent of space
        // Step 2: binning the space
        max_axis = cur_node->max_extent_axis(aabb_infos, bins, idx_bins);

        // Step 3: forward-backward linear sweep for heuristic calculation
        std::array<int, num_bins> prim_cnts;
        std::array<Ty, num_bins> fwd_areas, bwd_areas;
        for (int i = 0; i < num_bins; i++) {
            fwd_bound   += idx_bins[i].bound;
            prim_cnts[i] = idx_bins[i].prim_cnt;
            fwd_areas[i] = fwd_bound.area();
            if (i > 0) {
                bwd_bound += idx_bins[num_bins - i].bound;
                bwd_areas[num_bins - 1 - i] = bwd_bound.area();
            }
        }
        std::partial_sum(prim_cnts.begin(), prim_cnts.end(), prim_cnts.begin());

        // Step 4: use the calculated area to computed the segment boundary
        for (int i = 0; i < num_bins - 1; i++) {
            Ty cost = traverse_cost + node_inv_area * 
                (Ty(prim_cnts[i]) * fwd_areas[i] + (node_prim_cnt - (prim_cnts[i])) * bwd_areas[i]);
            if (cost < min_cost) {
                min_cost = cost;
                seg_bin_idx = i;
            }
        }
        fwd_bound.clear();
        bwd_bound.clear();
        for (int i = 0; i <= seg_bin_idx; i++)       // calculate child node bound
            fwd_bound += idx_bins[i].bound;
        for (int i = num_bins - 1; i > seg_bin_idx; i--)
            bwd_bound += idx_bins[i].bound;
        child_prim_cnt = prim_cnts[seg_bin_idx];
        split_pos      = bins[seg_bin_idx];
    }
   
    if (child_prim_cnt > 0 && min_cost < node_prim_cnt) {             // cost of splitting is less than making this node a leaf node
        // Step 5: split the node and initialize the children
        cur_node->split_leaf_node(aabb_infos, split_pos, max_axis, std::move(fwd_bound), std::move(bwd_bound));
        
        // Step 6: start recursive splitting for the children
        recursive_kd_tree_SAH(cur_node->lchild(), aabb_infos, depth + 1);
        recursive_kd_tree_SAH(cur_node->rchild(), aabb_infos, depth + 1);
    } else {
        // This is a leaf node, yet this is the only way that a leaf node contains more than one primitive
        cur_node->split_axis = NONE;
    }
}

template <typename Ty>
KDPrimNode<Ty>* kd_tree_root_start(const py::array_t<Ty>& world_min, const py::array_t<Ty>& world_max, int& node_num, std::vector<BVHInfo<Ty>>& bvh_infos) {
    // Build BVH tree root node and start recursive tree construction
    std::vector<int> index_vec;
    std::iota(index_vec.begin(), index_vec.end(), 0);
    KDPrimNode<Ty>* root_node = new KDPrimNode<Ty>(AABB<Ty>(), std::move(index_vec));
    // All the indices in the leaf nodes will be linearized as well
    Eigen::Vector3<Ty> &bound_min = root_node->aabb.mini, &bound_max = root_node->aabb.maxi;
    const Ty* const min_ptr = world_min.data(0), * const max_ptr = world_max.data(0);
    for (int i = 0; i < 3; i++) {
        bound_min(i) = min_ptr[i];
        bound_max(i) = max_ptr[i];
    }
    recursive_kd_tree_SAH(root_node, bvh_infos);
    return root_node;
}

// This is the final function call for `bvh_build`
template <typename Ty>
int recursive_linearize(KDPrimNode<Ty>* cur_node, std::vector<LinearKdNode<Ty>>& lin_nodes, std::vector<int>& indices, int& num_nodes) {
    lin_nodes.emplace_back(cur_node);
    int old_num = num_nodes;
    num_nodes ++;
    if (cur_node->is_leaf()) {
        lin_nodes.back().r_offset   = 0;
        lin_nodes.back().all_offset = 1;
        lin_nodes.back().idx_base = static_cast<int>(indices.size());
        lin_nodes.back().idx_num  = cur_node->num_elems;
        for (int idx: cur_node->sub_idxs)
            indices.push_back(idx);
    } else {
        lin_nodes.back().idx_base = 0;
        lin_nodes.back().idx_num  = 0;
        if (cur_node->lchild() != nullptr) {
            recursive_linearize(cur_node->lchild(), lin_nodes, std::vector<int>& indices, num_nodes);
            lin_nodes[old_num].r_offset = num_nodes - old_num;
        }
        if (cur_node->rchild() != nullptr) {
            recursive_linearize(cur_node->rchild(), lin_nodes, std::vector<int>& indices, num_nodes);
            lin_nodes[old_num].all_offset = num_nodes - old_num;
        }
    }
}

template<typename Ty>
SplitAxis KDPrimNode<Ty>::max_extent_axis(
    const std::vector<AABB<Ty>>& aabbs, 
    std::vector<Ty>& bins, std::vector<AxisBins<Ty>>& idx_bins
) const {
    Vec3 min_ctr, max_ctr;
    min_ctr.setConstant(1e9);
    max_ctr.setConstant(-1e9);
    for (auto index: *(this->sub_idxs)) {
        min_ctr = min_ctr.cwiseMin(aabbs[index].mini);
        max_ctr = min_ctr.cwiseMin(aabbs[index].maxi);
    }
    Vec3 diff = max_ctr - min_ctr;
    Ty max_diff = diff(0);
    int split_axis = 0;
    for (int i = 1; i < 3; i++) {
        if (diff(i) > max_diff) {
            max_diff = diff(i);
            split_axis = i;
        }
    }
    bins.resize(num_bins);
    Ty min_r = min_ctr(split_axis) - 0.001, interval = (max_diff + 0.002) / Ty(num_bins);
    for (size_t i = 0; i < bins.size(); i++) 
        bins[i] = min_r + interval * static_cast<Ty>(i);
    for (auto index: *(this->sub_idxs)) {
        size_t bin_idxs = std::lower_bound(bins.begin(), bins.end(), aabb_infos[index].axis_centroid(split_axis)) - bins.begin();
        idx_bins[bin_idxs].push(aabb_infos[index]);
    }
    return SplitAxis(split_axis);
}

template<typename Ty>
void KDPrimNode<Ty>::split_leaf_node(
    const std::vector<AABB<Ty>>& aabbs, Ty split_pos, SplitAxis split_axis,
    AABB<Ty>&& lchild_aabb, AABB<Ty>&& rchild_aabb
) {
    this->split_pos  = split_pos; 
    this->split_axis = split_axis; 

    std::vector<int> lcontainer, rcontainer;
    for (int idx: *sub_idxs) {
        if (aabbs[idx].axis_centroid(split_axis) < this->split_pos)
            lcontainer.push_back(idx);
        else
            rcontainer.push_back(idx);
    }
    sub_idxs.reset(nullptr);
    add_child(std::move(lchild_aabb), std::move(lcontainer), true);
    add_child(std::move(rchild_aabb), std::move(rcontainer), false);
}

template<typename Ty>
size_t KDPrimNode<Ty>::get_size() const {
    size_t size = sizeof(*this);  // Size of the current instance
    if (is_leaf())
        size += sub_idxs->size() * sizeof(int);
    size += sizeof(_lchild);
    size += sizeof(_rchild);
    return size;
}

template class KDPrimNode<float>;               // 3D KD-tree (vector container)
template class KDPrimNode<double>;              // 3D KD-tree (vector container)

} // namespace scds