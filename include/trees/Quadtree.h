#pragma once
/**
 * Quad tree CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/
#include "Tree.h"
#include <vector>

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

protected:
    Point2<T> tl_pt;
    Point2<T> br_pt;
    Point2<T> tree_volumn;
};

}       // end namespace scds