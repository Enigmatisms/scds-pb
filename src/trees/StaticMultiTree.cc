#include <queue>
#include "trees/StaticMultiTree.h"

namespace scds {

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::search_nn_bf(const Pointx& pt, PointVec& nn, size_t k, T radius) const {
    auto distance_comp = \
    [](const auto& pr1, const auto& pr2) {
        return pr1.second < pr2.second;
    };
    auto radius2 = radius * radius;
    std::priority_queue<std::pair<size_t, T>, std::vector<std::pair<size_t, T>>, decltype(distance_comp)> max_heap(distance_comp);
    for (size_t pt_idx: root->get_indices()) {
        auto query_p = (*all_pts)[pt_idx];
        T distance2  = (query_p - pt).length2();
        if (distance2 > radius2) continue;
        
        max_heap.emplace(pt_idx, distance2);
        if (max_heap.size() > k) 
            max_heap.pop();
    }
    while (!max_heap.empty()) {
        size_t idx = max_heap.top().first;
        max_heap.pop();
        nn.push_back((*all_pts)[idx]);
    }
}

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::search_nn(const Pointx& pt, PointVec& nn, size_t k, T radius) const {
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

    auto search_tr = pt + float(radius), search_bl = pt - float(radius);
    T radius2 = radius * radius;
    std::vector<NodePtr> stack;
    stack.reserve(32);
    stack.push_back(root);

    auto distance_comp = \
    [](const auto& pr1, const auto& pr2) {
        return pr1.second < pr2.second;
    };
    std::priority_queue<std::pair<size_t, T>, std::vector<std::pair<size_t, T>>, decltype(distance_comp)> max_heap(distance_comp);
    
    while (!stack.empty()) {
        NodePtr top_node = stack.back();
        stack.pop_back();
        for (size_t i = 0; i < Nchild; i++) {
            auto child = top_node->get_child(i);       
            if (!child || !child->overlap_range(search_tr, search_bl)) continue;
            if (!child->is_leaf) {
                stack.push_back(child);
                continue;
            }
            for (size_t pt_idx: child->get_indices()) {
                auto query_p = (*all_pts)[pt_idx];      // assign constructor seems to be faulty
                T distance2  = (query_p - pt).length2();
                if (distance2 > radius2) continue;
                
                max_heap.emplace(pt_idx, distance2);
                if (max_heap.size() > k) 
                    max_heap.pop();
            }
        }
    }
    while (!max_heap.empty()) {
        size_t idx = max_heap.top().first;
        max_heap.pop();
        nn.push_back((*all_pts)[idx]);
    }
}

template<typename T, size_t Ndim, size_t Nchild>
void StaticMultiTree<T, Ndim, Nchild>::insert(const Pointx& pt) {
    size_t new_index = all_pts->size(), cur_depth = 0;
    all_pts->emplace_back(pt);
    auto ptr = root;
    do {
        ptr->insert(new_index);         // index will be inserted imediately
        // if we can (and must, since some condition is violated) built sub-trees, then:
        if (cur_depth < max_depth && ptr->num_points() >= node_max_point_num) {
            if (!ptr->is_leaf) {
                // current note is not leaf, decide which quandrant the point is in
                size_t child_id = which_child(ptr, pt);
                ptr = ptr->try_get_child(child_id);
                cur_depth ++;
                continue;
            }
            // parition the tree
            ptr->is_leaf = false;
            bool same_child = false;
            do {
                std::array<std::unordered_set<size_t>, Nchild> sub_sets;
                for (auto idx: ptr->get_indices()) {
                    // we should consider the worst case: that the N + 1 points recursively fail to be 
                    // partioned (since they are close to each other and always fall in the same quadrant)  
                    const Pointx& p = (*all_pts)[idx];
                    size_t child_id = which_child(ptr, p);
                    sub_sets[child_id].emplace(idx);
                }
                size_t next_child_id = 0;
                for (size_t i = 0; i < Nchild; i++) {
                    if (sub_sets.size() > node_max_point_num) {
                        same_child = true;
                        next_child_id  = i;
                    }
                }
                if (!same_child) {
                    // points are not in the same quadrant, the leaf node can be successfully partitioned
                    Pointx half_size = ptr->size / 2, offset;
                    for (size_t child_id = 0; child_id < Nchild; child_id++) {
                        if (sub_sets[child_id].empty()) continue;

                        auto offset = Node::get_child_offset(half_size, child_id);
                        ptr->add_child(ptr->center + offset, half_size, std::move(sub_sets[child_id]), child_id);
                    }
                    return;
                } else {
                    // update ptr to that same
                    if (cur_depth >= max_depth) {
                        // since max depth is reached, no child can be created, we are making the current node a leaf node
                        ptr->is_leaf = true;   
                        return;
                    }
                    ptr = ptr->try_get_child(next_child_id);
                    cur_depth ++;
                }
                // check whether all the points are in the same quadrant
            } while (same_child);
        } else {                    // no need to build sub-tree
            return;
        }
    } while (true);
}

template class StaticMultiTree<float, 2, 4>;
template class StaticMultiTree<float, 3, 8>;

}   // end namespace scds