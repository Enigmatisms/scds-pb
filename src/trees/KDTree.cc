#include <queue>
#include "utils/utils.h"
#include "utils/pybind_utils.h"
#include "trees/KDTree.h"

namespace scds {

template<typename T, size_t Ndim>
void KDTree<T, Ndim>::recursive_solve(const Pointx& pt, QueueType& queue, std::shared_ptr<Node> cur_node) const {
    if (!cur_node->is_leaf()) {
        auto res_child = cur_node->resident_child(pt);
        recursive_solve(pt, queue, res_child);
        T axial_diff = pt[cur_node->split_axis] - cur_node->split_pos;
        axial_diff *= axial_diff;
        if ((!queue.empty() && queue.top().second > axial_diff) || (static_cast<int>(queue.size()) < k && axial_diff < radius2)) {
            recursive_solve(pt, queue, cur_node->the_other(res_child.get()));
        }
    } else {        // leaf node: traverse all the points in the node
        for (int pt_idx: cur_node->get_indices()) {
            T distance2  = ((*all_pts)[pt_idx] - pt).length2();
            if (distance2 > radius2) continue;
            if (static_cast<int>(queue.size()) >= k) {
                if (distance2 < queue.top().second) {
                    queue.pop();
                    queue.emplace(pt_idx, distance2);
                }
            } else {
                queue.emplace(pt_idx, distance2);
            }
        }
    }
}

template<typename T, size_t Ndim>
void KDTree<T, Ndim>::search_nn(const Pointx& pt, PointVec& nn, int k, T radius) const {
    ProfilePhase _(Prof::KDTreeSearchNN);
    if (!k) {
        printf("Warning: specified nearest neight num is 0. Set 1 by default.\n");
        k = 1;
    }
    this->k = k;
    this->radius = radius;
    this->radius2 = radius * radius;

    QueueType max_heap(DistanceComp<T>{});

    recursive_solve(pt, max_heap, root);

    nn.reserve(k);
    while (!max_heap.empty()) {
        int idx = max_heap.top().first;
        max_heap.pop();
        nn.push_back((*all_pts)[idx]);
    }
}

template<typename T, size_t Ndim>
void KDTree<T, Ndim>::insert(const Pointx& pt) {
    ProfilePhase _(Prof::StaticMultiTreeInsert);
    int new_index = static_cast<int>(all_pts->size()), cur_depth = 0;
    all_pts->push_back(pt);
    auto ptr = root;

    #ifdef SMT_MEMORY_PROFILE
        size_t pt_bytes = sizeof(Pointx);
        pointBytes += pt_bytes;
    #endif //SMT_MEMORY_PROFILE
    // This incremental tree construction will be very inefficient (for static multi-tree, it is different)
    do {
        ptr->num_points++;
        if (cur_depth < max_depth && ptr->num_points > node_max_point_num) {
            tree_depth = std::max(tree_depth, cur_depth ++) + 1;
            if (!ptr->is_leaf()) {
                // find the child leaf the new point should reside in
                ptr = ptr->resident_child(pt);
                continue;
            }
            // When the current node should be splitted (maximum point riched)
            ptr->insert(new_index);
            ptr->split_leaf_node(*all_pts);
            return;
        } else {
            ptr->insert(new_index);
            return;
        }
    } while (true);
    
}

template<typename T, size_t Ndim>
void KDTree<T, Ndim>::build_tree_py(const pybind11::array_t<T>& pts) {
    const size_t num_pts = pts.shape()[0];
    std::vector<Pointx> points;
    points.reserve(num_pts);
    const T* ptr = pts.data();
    for (size_t i = 0; i < num_pts; i++, ptr += Ndim)
        points.emplace_back(Pointx::from_pointer(ptr));
    build_tree(std::move(points));
}

template<typename T, size_t Ndim>
void KDTree<T, Ndim>::insert_py(const pybind11::array_t<T>& pt) {
    ProfilePhase _(Prof::StaticMultiTreeInsertPy);
    SMT_ARRAY_DTYPE_CHECK(pt, T);
    bool is_single = pyArrayShapeCheck(pt, Ndim);
    if (is_single) {
        insert(Pointx::from_pointer(pt.data()));
    } else {
        const T* ptr = pt.data();
        const size_t num_pts = pt.shape()[0];
        for (size_t i = 0; i < num_pts; i++, ptr += Ndim)
            insert(Pointx::from_pointer(ptr));
    }
}

template<typename T, size_t Ndim>
pybind11::array_t<T> KDTree<T, Ndim>::search_nn_py(const pybind11::array_t<T>& pt, int k, T radius) const {
    ProfilePhase _(Prof::StaticMultiTreeSearchNNPy);

    SMT_ARRAY_DTYPE_CHECK(pt, T);
    bool is_single = pyArrayShapeCheck(pt, Ndim);
    if (!is_single)
        SCDS_RUNTIME_ERROR("`search` function can only search one point at a time.");

    auto to_search = Pointx::from_pointer(pt.data());
    std::vector<Pointx> nn;
    search_nn(to_search, nn, k, radius);

    pybind11::array_t<T> result = create_array_2d<T>(nn.size(), Ndim);
    T* ptr = result.mutable_data();
    for (size_t i = 0; i < nn.size(); i++, ptr+=Ndim)
        memcpy(ptr, nn[i].const_data(), sizeof(T) * Ndim);
    return result;
}

template<typename T, size_t Ndim>
pybind11::tuple KDTree<T, Ndim>::get_tree_structure() const {
    static_assert(Ndim < 4, "Visualizing 4+ dimension is pointless.");
    using NodePtr = std::shared_ptr<Node>;
    std::vector<NodePtr> stack;
    stack.reserve(32);
    stack.push_back(root);
    constexpr size_t Ndim2 = Ndim << 1;

    std::vector<Point4<T>> non_leaf_nodes, leaf_nodes;
    non_leaf_nodes.reserve(root->num_points << 1);
    leaf_nodes.reserve(root->num_points);

    while (!stack.empty()) {
        NodePtr top_node = stack.back();
        Point<T, Ndim2> node_range;
        for (size_t i = 0; i < Ndim; i++) {
            node_range[i]        = top_node->center[i];
            node_range[i + Ndim] = top_node->half_size[i];
        }
        if (top_node->is_leaf())
            leaf_nodes.push_back(node_range);
        else
            non_leaf_nodes.push_back(node_range);
        stack.pop_back();
        auto lchild = top_node->lchild(),       
             rchild = top_node->rchild();       
        if (lchild) stack.push_back(lchild);
        if (rchild) stack.push_back(rchild);
    }
    pybind11::array_t<T> non_leaf_nodes_py = create_array_2d<T>(non_leaf_nodes.size(), Ndim2);
    pybind11::array_t<T> leaf_nodes_py     = create_array_2d<T>(leaf_nodes.size(), Ndim2);
    T* ptr1 = non_leaf_nodes_py.mutable_data(), *ptr2 = leaf_nodes_py.mutable_data();
    for (size_t i = 0; i < non_leaf_nodes.size(); i++, ptr1+=Ndim2)
        memcpy(ptr1, non_leaf_nodes[i].const_data(), sizeof(T) * Ndim2);
    for (size_t i = 0; i < leaf_nodes.size(); i++, ptr2+=Ndim2)
        memcpy(ptr2, leaf_nodes[i].const_data(), sizeof(T) * Ndim2);
    return pybind11::make_tuple(non_leaf_nodes_py, leaf_nodes_py);
}


template class KDTree<float, 2>;
template class KDTree<float, 3>;
template class KDTree<double, 2>;
template class KDTree<double, 3>;

}   // end namespace scds