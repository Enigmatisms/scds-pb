#include <queue>
#include "utils/utils.h"
#include "utils/pybind_utils.h"
#include "trees/StaticMultiTree.h"

namespace scds {

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::search_nn_bf(const Pointx& pt, PointVec& nn, size_t k, T radius) const {
    ProfilePhase _(Prof::StaticMultiTreeSearchNNBF);
    auto distance_comp = \
    [](const auto& pr1, const auto& pr2) {
        return pr1.second < pr2.second;
    };
    auto radius2 = radius * radius;
    std::priority_queue<std::pair<int, T>, std::vector<std::pair<int, T>>, decltype(distance_comp)> max_heap(distance_comp);
    for (int pt_idx = 0; pt_idx < root->num_points; pt_idx++) {
        auto query_p = (*all_pts)[pt_idx];
        T distance2  = (query_p - pt).length2();
        if (distance2 > radius2) continue;
        
        if (max_heap.size() >= k) {
            if (distance2 < max_heap.top().second) {
                max_heap.pop();
                max_heap.emplace(pt_idx, distance2);
            }
        } else {
            max_heap.emplace(pt_idx, distance2);
        }
    }
    nn.reserve(k);
    while (!max_heap.empty()) {
        int idx = max_heap.top().first;
        max_heap.pop();
        nn.push_back((*all_pts)[idx]);
    }
}

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::search_nn(const Pointx& pt, PointVec& nn, size_t k, T radius) const {
    ProfilePhase _(Prof::StaticMultiTreeSearchNN);
    using NodePtr = std::shared_ptr<Node>;

    if (radius < 1e-5) {
        printf("Warning: Search range not specified, radius = %f\n", radius);
        printf("You should pass in a reasonable search radius.\n");
        return;
    }
    if (!k) {
        printf("Warning: specified nearest neight num is 0. Set 1 by default.\n");
        k = 1;
    }

    auto search_tr = pt + static_cast<float>(radius), search_bl = pt - static_cast<float>(radius);
    T radius2 = radius * radius;
    std::vector<NodePtr> stack;
    stack.reserve(32);
    stack.push_back(root);

    auto distance_comp = \
    [](const auto& pr1, const auto& pr2) {
        return pr1.second < pr2.second;
    };
    std::priority_queue<std::pair<int, T>, std::vector<std::pair<int, T>>, decltype(distance_comp)> max_heap(distance_comp);
    while (!stack.empty()) {
        NodePtr top_node = stack.back();
        stack.pop_back();
        for (size_t i = 0; i < Nchild; i++) {
            auto child = top_node->get_child(i);       
            if (!child || !child->overlap_range(search_tr, search_bl)) continue;
            if (!child->is_leaf()) {
                stack.push_back(child);
                continue;
            }
            for (int pt_idx: child->get_indices()) {
                T distance2  = ((*all_pts)[pt_idx] - pt).length2();
                if (distance2 > radius2) continue;
                if (max_heap.size() >= k) {
                    if (distance2 < max_heap.top().second) {
                        max_heap.pop();
                        max_heap.emplace(pt_idx, distance2);
                    }
                } else {
                    max_heap.emplace(pt_idx, distance2);
                }
            }
        }
    }
    nn.reserve(k);
    while (!max_heap.empty()) {
        int idx = max_heap.top().first;
        max_heap.pop();
        nn.push_back((*all_pts)[idx]);
    }
}

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::insert(const Pointx& pt) {
    ProfilePhase _(Prof::StaticMultiTreeInsert);
    int new_index = static_cast<int>(all_pts->size()), cur_depth = 0;
    all_pts->push_back(pt);
    auto ptr = root;

    #ifdef SMT_MEMORY_PROFILE
        size_t pt_bytes = sizeof(Pointx);
        pointBytes += pt_bytes;
    #endif //SMT_MEMORY_PROFILE

    do {
        ptr->num_points++;
        // if we can (and must, since some condition is violated) built sub-trees, then:
        if (cur_depth < max_depth && ptr->num_points > node_max_point_num) {
            if (!ptr->is_leaf()) {
                // current note is not leaf, decide which quandrant the point is in
                size_t child_id = which_child(ptr, pt);
                ptr = ptr->try_get_child(child_id);
                tree_depth = std::max(tree_depth, cur_depth ++) + 1;
                continue;
            }
            ptr->insert(new_index);
            bool same_child = false;
            do {
                std::array<std::vector<int>, Nchild> sub_sets;
                for (auto idx: ptr->get_indices()) {
                    // we should consider the worst case: that the N + 1 points recursively fail to be 
                    // partioned (since they are close to each other and always fall in the same quadrant)  
                    const Pointx& p = (*all_pts)[idx];
                    size_t child_id = which_child(ptr, p);
                    sub_sets[child_id].emplace_back(idx);
                }

                same_child = false;
                size_t next_child_id = 0;
                for (size_t i = 0; i < Nchild; i++) {
                    if (static_cast<int>(sub_sets[i].size()) > node_max_point_num) {
                        same_child = true;
                        next_child_id  = i;
                    }
                }
                
                if (!same_child) {
                    #ifdef TREE_NODE_MEMORY_PROFILE
                        leafNodes --;
                    #endif //TREE_NODE_MEMORY_PROFILE
                    ptr->set_non_leaf();
                    // points are not in the same quadrant, the leaf node can be successfully partitioned
                    Pointx half_size = ptr->size / 2, offset;
                    for (size_t child_id = 0; child_id < Nchild; child_id++) {
                        if (sub_sets[child_id].empty()) continue;

                        auto offset = Node::get_child_offset(half_size, child_id);
                        ptr->add_child(ptr->center + offset, half_size, std::move(sub_sets[child_id]), child_id);
                    }
                    return;
                } else {
                    if (cur_depth >= max_depth) {
                        // since max depth is reached, no child can be created, we are making the current node a leaf node
                        return;
                    }
                    #ifdef TREE_NODE_MEMORY_PROFILE
                        leafNodes --;
                    #endif //TREE_NODE_MEMORY_PROFILE
                    ptr->set_non_leaf();   
                    ptr = ptr->try_get_child(next_child_id);
                    ptr->overwrite_sub_idxs(std::move(sub_sets[next_child_id]));
                    tree_depth = std::max(tree_depth, cur_depth ++) + 1;
                }
                // check whether all the points are in the same quadrant
            } while (same_child);
        } else {                    // no need to build sub-tree
            ptr->insert(new_index);
            return;
        }
    } while (true);
}

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::insert_py(const pybind11::array_t<T>& pt) {
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

template<typename T, size_t Ndim, size_t Nchild>
pybind11::array_t<T> StaticMultiTree<T, Ndim, Nchild>::search_nn_py(const pybind11::array_t<T>& pt, int k, T radius) const {
    ProfilePhase _(Prof::StaticMultiTreeSearchNNPy);

    SMT_ARRAY_DTYPE_CHECK(pt, T);
    bool is_single = pyArrayShapeCheck(pt, Ndim);
    if (!is_single)
        SCDS_RUNTIME_ERROR("`search` function can only search one point at a time.");
    auto to_search = Pointx::from_pointer(pt.data());
    std::vector<Pointx> nn;
    search_nn(to_search, nn, k, radius);
    // TODO: this can be optimized - copying point data from PointVec to array_t
    pybind11::array_t<T> result = create_array_2d<T>(nn.size(), Ndim);
    T* ptr = result.mutable_data();
    for (size_t i = 0; i < nn.size(); i++, ptr+=Ndim)
        memcpy(ptr, nn[i].const_data(), sizeof(T) * Ndim);
    return result;
}

template<typename T, size_t Ndim, size_t Nchild>
pybind11::array_t<T> StaticMultiTree<T, Ndim, Nchild>::search_nn_bf_py(const pybind11::array_t<T>& pt, int k, T radius) const {
    ProfilePhase _(Prof::StaticMultiTreeSearchNNBFPy);

    SMT_ARRAY_DTYPE_CHECK(pt, T);
    bool is_single = pyArrayShapeCheck(pt, Ndim);
    if (!is_single)
        SCDS_RUNTIME_ERROR("`search` function can only search one point at a time.");
    auto to_search = Pointx::from_pointer(pt.data());
    std::vector<Pointx> nn;
    search_nn_bf(to_search, nn, k, radius);

    // TODO: this can be optimized - copying point data from PointVec to array_t
    pybind11::array_t<T> result = create_array_2d<T>(nn.size(), Ndim);
    T* ptr = result.mutable_data();
    for (size_t i = 0; i < nn.size(); i++, ptr+=Ndim)
        memcpy(ptr, nn[i].const_data(), sizeof(T) * Ndim);
    
    return result;
}

template<typename T, size_t Ndim, size_t Nchild>
pybind11::tuple StaticMultiTree<T, Ndim, Nchild>::get_tree_structure() const {
    static_assert(Ndim < 4, "Visualizing 4+ dimension is point less.");
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
            node_range[i + Ndim] = top_node->size[i];
        }
        if (top_node->is_leaf())
            leaf_nodes.push_back(node_range);
        else
            non_leaf_nodes.push_back(node_range);
        stack.pop_back();
        for (size_t i = 0; i < Nchild; i++) {
            auto child = top_node->get_child(i);       
            if (!child) continue;
            stack.push_back(child);
        }
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


template class StaticMultiTree<float, 2, 4>;
template class StaticMultiTree<float, 3, 8>;
template class StaticMultiTree<double, 2, 4>;
template class StaticMultiTree<double, 3, 8>;

}   // end namespace scds