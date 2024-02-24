#include "prims/KDPrims.h"

namespace scds {

constexpr int num_bins = 1;
constexpr int split_node_prim_num = 8;
constexpr float traverse_cost = 0.1;

template <typename Ty>
int recursive_kd_tree_SAH(KDPrimNode<Ty>* const cur_node, std::vector<AABB<Ty>>& aabb_infos) {
    AABB fwd_bound, bwd_bound;
    int seg_idx = 0, child_prim_cnt = 0;                // this index is used for indexing variable `bins`
    const int prim_num = cur_node->prim_num, base = cur_node->base, max_pos = base + prim_num;
    Ty min_cost = 5e9, node_prim_cnt = Ty(prim_num), node_inv_area = 1. / cur_node->bound.area();

    // Step 1: decide the axis that expands the maximum extent of space
    std::vector<Ty> bins;        // bins: from (start_pos + interval) to end_pos
    // FIXME: the bin bound calculation
    SplitAxis max_axis = cur_node->max_extent_axis(aabb_infos, bins);
    if (cur_node->prim_num > split_node_prim_num) {   // SAH

        // Step 2: binning the space
        std::array<AxisBins, num_bins> idx_bins;
        for (int i = cur_node->base; i < max_pos; i++) {
            size_t index = std::lower_bound(bins.begin(), bins.end(), aabb_infos[i].axis_centroid(max_axis)) - bins.begin();
            idx_bins[index].push(aabb_infos[i]);
        }

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
        int seg_bin_idx = 0;
        for (int i = 0; i < num_bins - 1; i++) {
            Ty cost = traverse_cost + node_inv_area * 
                (Ty(prim_cnts[i]) * fwd_areas[i] + (node_prim_cnt - (prim_cnts[i])) * bwd_areas[i]);
            if (cost < min_cost) {
                min_cost = cost;
                seg_bin_idx = i;
            }
        }
        // Step 5: reordering the BVH info in the vector to make the segment contiguous (partition around pivot)
        // FIXME: also, impose the depth constraints
        if (min_cost < node_prim_cnt) {
            // FIXME: here we break the sub_idxs of the current node, and create two child node?
            // FIXME: we don't need to keep track of the offsets? which can be computed only
            // FIXME: This is not partition here, since one primitive can be observed in bith lchild and rchild
            std::partition(aabb_infos.begin() + base, aabb_infos.begin() + max_pos,
                [pivot = bins[seg_bin_idx], dim = max_axis](const BVHInfo& bvh) {
                    return bvh.centroid[dim] < pivot;
            });
            // FIXME: the following might be preserved
            child_prim_cnt = prim_cnts[seg_bin_idx];
        }
        fwd_bound.clear();
        bwd_bound.clear();
        for (int i = 0; i <= seg_bin_idx; i++)       // calculate child node bound
            fwd_bound += idx_bins[i].bound;
        for (int i = num_bins - 1; i > seg_bin_idx; i--)
            bwd_bound += idx_bins[i].bound;
    } else {                                    // equal primitive number
        seg_idx = (base + max_pos) >> 1;
        // Step 5: reordering the BVH info in the vector to make the segment contiguous (keep around half of the bvh in lchild)
        std::nth_element(aabb_infos.begin() + base, aabb_infos.begin() + seg_idx, aabb_infos.begin() + max_pos,
            [dim = max_axis] (const BVHInfo& bvh1, const BVHInfo& bvh2) {
                return bvh1.centroid[dim] < bvh2.centroid[dim];
            }
        );
        for (int i = base; i < seg_idx; i++)    // calculate child node bound
            fwd_bound += aabb_infos[i].bound;
        for (int i = seg_idx; i < max_pos; i++)
            bwd_bound += aabb_infos[i].bound;
        child_prim_cnt = seg_idx - base;        // bvh[seg_idx] will be in rchild
        Ty split_cost = traverse_cost + node_inv_area * 
                (fwd_bound.area() * child_prim_cnt + bwd_bound.area() * (node_prim_cnt - child_prim_cnt));
        if (split_cost >= node_prim_cnt)
            child_prim_cnt = 0;
    }
    if (child_prim_cnt > 0) {             // cost of splitting is less than making this node a leaf node
        // this should be fixed
        // Step 5: split the node and initialize the children
        cur_node->lchild = new BVHNode(base, child_prim_cnt);
        cur_node->rchild = new BVHNode(base + child_prim_cnt, prim_num - child_prim_cnt);

        cur_node->lchild->bound = fwd_bound;
        cur_node->rchild->bound = bwd_bound;
        cur_node->axis = max_axis;
        // Step 7: start recursive splitting for the children
        int node_num = 1;
        if (cur_node->lchild->prim_num > max_node_prim)
            node_num += recursive_bvh_SAH(cur_node->lchild, aabb_infos);
        else node_num ++;
        if (cur_node->rchild->prim_num > max_node_prim)
            node_num += recursive_bvh_SAH(cur_node->rchild, aabb_infos);
        else node_num ++;
        return node_num;
    } else {
        // This is a leaf node, yet this is the only way that a leaf node contains more than one primitive
        cur_node->axis = NONE;
        return 1;
    }
}

template<typename Ty>
SplitAxis KDPrimNode<Ty>::max_extent_axis(const std::vector<AABB<Ty>>& pts) const {
    return SplitAxis::NONE;
}

template<typename Ty>
void KDPrimNode<Ty>::split_leaf_node(const std::vector<Eigen::Vector3<Ty>>& pts) {
    // get split axis
    SplitAxis axis   = max_extent_axis(pts);
    this->split_axis = axis; 
    // get split pos using nth_element
    auto comp_op = [axis, &pts](int index_1, int index_2) {return pts[index_1][axis] < pts[index_2][axis];};

    int half_num = static_cast<int>(sub_idxs->size() - 1) >> 1;
    std::nth_element(sub_idxs->begin(), sub_idxs->begin() + half_num, sub_idxs->end(), comp_op);
    int index_1 = *(sub_idxs->begin() + half_num);
    std::nth_element(sub_idxs->begin(), sub_idxs->begin() + half_num + 1, sub_idxs->end(), comp_op);
    int index_2 = *(sub_idxs->begin() + half_num + 1);

    this->split_pos = (pts[index_1][axis] + pts[index_2][axis]) * static_cast<Ty>(0.5);

    std::vector<int> lcontainer, rcontainer;
    for (int idx: *sub_idxs) {
        if (pts[idx][axis] < this->split_pos)
            lcontainer.push_back(idx);
        else
            rcontainer.push_back(idx);
    }
    sub_idxs.reset(nullptr);
    Ty lsize = (half_size[axis] - center[axis] + split_pos) / 2, rsize = half_size[axis] - lsize,
       lpos  = split_pos - lsize, rpos = split_pos + rsize;
    Pointx l_ctr = this->center, r_ctr = this->center,
           l_size = this->half_size, r_size = this->half_size;
    l_ctr[axis]  = lpos;
    l_size[axis] = lsize;
    r_ctr[axis]  = rpos;
    r_size[axis] = rsize;
    add_child(std::move(l_ctr), std::move(l_size), std::move(lcontainer), true);
    add_child(std::move(r_ctr), std::move(r_size), std::move(rcontainer), false);
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