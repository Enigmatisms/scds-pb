#pragma once
/**
 * Quad tree CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/
#include "tree.h"

namespace scds {

/**
 * Static quadtree for static scenes. Since quad trees work in
 * 2D spaces, it will be much easier to visualize and debug
 * 
 * Note that quadtree and octree are not frequently used to perform NN search
 * Quadtree and octree are used to subdivide the space to store some specific information
 * like occupancy grid (in SLAM) and radiance (in rendering - path guiding)
*/
class StaticQuadTree: public TreeBase {
public:
    StaticQuadTree();
};

}       // end namespace scds