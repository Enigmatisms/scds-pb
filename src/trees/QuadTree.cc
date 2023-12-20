#include "trees/Quadtree.h"

namespace scds {

template<typename T>
template<typename PointType>
void StaticQuadTree<T>::insert(PointType&& pt) {
    ASSERT_POINT_TYPE(PointType);
    // Loop (until the break condition is met)
    // step 1: whether the current node is leaf / empty. If is leaf: to step 2, else to step 3
    // step 2: add to the current leaf node, then check the tree partition condition
    // if the tree should be further partitioned, for every pt in the subtree, decide which quadrant each note is in
    size_t new_index = all_pts.size();
    all_pts.push_back(pt);
    QdNode<T>* ptr = root.get();
    size_t cur_depth = 0;
    do {
        ptr->insert(new_index);         // index will be inserted imediately
        // decide quadrant
        // Depth can be extended, also number of points reaches the maximum
        if (cur_depth < max_depth && ptr->num_points() >= node_max_point_num) {
            if (ptr->is_leaf) {     // parition the tree
                ptr->is_leaf = false;
                bool same_quadrant = false;
                do {
                    std::array<std::unordered_set<size_t>, Ndim> sub_sets;
                    for (auto idx: ptr->sub_idxs) {
                        // we should consider the worst case: that the N + 1 points recursively fail to be 
                        // partioned (since they are close to each other and always fall in the same quadrant)  
                        const Point2<T>& p = all_pts[idx];
                        size_t quad_id = which_quadrant(ptr, p);
                        sub_sets[quad_id].emplace(idx);
                    }
                    size_t next_quad_id = 0;
                    for (size_t i = 0; i < Ndim; i++) {
                        if (sub_sets.size() > node_max_point_num) {
                            same_quadrant = true;
                            next_quad_id  = i;
                        }
                    }
                    if (!same_quadrant) {
                        // points are not in the same quadrant, the leaf node can be successfully partitioned
                        Point2<T> half_size = ptr->size / 2, offset;
                        for (size_t quad_id = 0; quad_id < Ndim; quad_id++) {
                            if (sub_sets[i].empty()) continue;

                            auto offset = QdNode<T>::quadrant_offset(half_size, quad_id);
                            ptr->add_child(ptr->center + offset, half_size, std::move(sub_sets[quad_id]), quad_id);
                        }
                        return;
                    } else {
                        // update ptr to that same
                        if (cur_depth < max_depth) {
                            ptr = ptr->get_child(next_quad_id);
                            cur_depth ++;
                        } else {
                            // since max depth is reached, no child can be created, we are making the current node a leaf node
                            ptr->is_leaf = true;   
                            return;
                        }
                    }
                    // check whether all the points are in the same quadrant
                } while (same_quadrant)
            } else {                // go to the sub-tree
                size_t quad_id = which_quadrant(ptr, pt);
                ptr = ptr->get_child(quad_id);
                cur_depth ++;
            }
        } else {
            return;
        }
    } while (true);
}

}   // end namespace scds