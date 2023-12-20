#pragma once
/**
 * Quad tree CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/
#include "Tree.h"
#include "utils/TreeNode.h"

namespace scds {

/**
 * Static quadtree for static scenes. Since quad trees work in
 * 2D spaces, it will be much easier to visualize and debug
 * 
 * Note that quadtree and octree are not frequently used to perform NN search
 * Quadtree and octree are used to subdivide the space to store some specific information
 * like occupancy grid (in SLAM) and radiance (in rendering - path guiding)
*/

template <typename T>
class StaticQuadTree: public TreeBase<T> {
public:
    template <typename U = T>
    using PointVec = std::vector<Point2<U>>;
    
    // lazy build (incremental)
    StaticQuadTree(const Point2<T>& tl, const Point2<T>& br, int max_depth = 0);
    // full build
    StaticQuadTree(const PointVec<T>& points, int max_depth = 0);
public:
    template <typename PointType>
    virtual void insert(PointType&& pt);
protected:
    /**
     * @brief decide which quadrant the node is in (via bit operation)
     *    y
     * 01 | 11 
     * ---.--- x
     * 00 | 10
    */
    template <typename PointType>
    static size_t which_quadrant(QdNode<T>* node, PointType&& pt) {
        ASSERT_POINT_TYPE(PointType);
        auto judge = ((pt - node->center) > 0).to_u64();
        return (judge[0] << 1) + judge[1];
    }
protected:
    Point2<T> tl_pt;
    Point2<T> br_pt;
    Point2<T> tree_volumn;
    PointVec<T> all_pts;
    std::shared_ptr<QdNode<T>> root;

private:
    // maximum depth of the tree (1st priority)
    size_t max_depth;
    // maximum number of point in a node (2nd priority)
    size_t node_max_point_num;

    static constexpr size_t MAX_DEPTH = 32;                         // physical barrier
};

}       // end namespace scds