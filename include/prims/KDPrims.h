#pragma once
/**
 * KD-tree for primitives. This will be migrated to my AdaPT renderer
 * @author: Qianyue He
 * @date:   2024-2-22
*/

#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <queue>
#include "utils/BinTreeNode.h"
#include "utils/stats.h"

namespace scds {

#define KDPT_MEMORY_PROFILE

#ifdef KDPT_MEMORY_PROFILE

STAT_MEMORY_COUNTER("KDPT/Total points", pointBytes);

#endif //KDPT_MEMORY_PROFILE

/**
 * @brief KD-tree efficient nearest neighbor earch structure
*/
template <typename T>
class KDPrimTree {
public:
    /**
     * main logic walk-through
     * 
     * special note:
     * (1) the tree construction will be a recursive process: primitives are indexed with indices
     * for the lchild and rchild of a node, should we re-order the primitives to accelerate querying (cache)?
     * I don't think so, since one primitive can exist in different nodes, this re-ordering might 
     * be tricky? Maybe it is only feasible for BVH.
     * (2) the splitting plane is chosen to be one of the planes of the bounding box
     * (3) the tree nodes will be linearized in a pre-order way
     * 
     * tree construction logic:
     * 1. for the current node, find the maximum extent axis (of bounding boxes), denoted by T
     * 2. along that axis, find the splitting plane (among the defining planes of the AABBs)
     * via SAH
     * 3. since AABBs can exist in both the lchild and rchild of a node, std::partition seems to be... meaningless
     * So what we do here is simply recording the indices to the AABBs to the child nodes (and note that, only leaf
     * node keep track of these indices. The indices field (unique_ptr) will be nullified if the current node is not leaf
     * 4. recursively split the leaf node until stopping criteria are met
     * 
     * ray intersection logic:
     * (1) Note that, before this is used in Taichi end, node linearization will be performed (pre-order)
     * We should be clear whether pre-order linearization or in-order linearization is better (since PBR-book performs
     * in-order linearization). Here we assume the linearization is pre-order one, therefore (root)->(l-tree)->(r-tree).
     * 
     * (2) This runs into the same pit-fall: in Taichi it will be un-easy to do optimal ray-intersection test, since
     * stack (todo-list) will not be easy to implement, unless we are using tile-based rendering (no-preview): pre-allocated
     * stack (with fixed upper limit of stack capacity) can be reused between different tiles. Otherwise, we should find a way
     * to allocate a local / variable array
     * 
     * Other parts are yet to be cleared
    */
};

}       // end namespace scds