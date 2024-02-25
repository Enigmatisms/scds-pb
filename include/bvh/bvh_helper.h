/**
 * @file bvh_helper.h
 * @author Qianyue He
 * @brief BVH utilities
 * @copyright Copyright (c) 2023
 */

#pragma once
#include <algorithm>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "utils/utils.h"


namespace scds {

namespace py = pybind11;
template <typename Ty>
struct AABB {
    Eigen::Vector3<Ty> mini;
    Eigen::Vector3<Ty> maxi;

    AABB() {
        mini.fill(1e4);
        maxi.fill(-1e4);
    }
    /** Review needed: is this good to implement two overload? */
    template <typename Ptype1, typename Ptype2>
    AABB(Ptype1&& mini, Ptype2&& maxi): mini(std::forward<Ptype1>(mini)), maxi(std::forward<Ptype2>(maxi)) {}
    AABB(const Eigen::Matrix3<Ty>& primitive, bool is_sphere = false) {
        if (is_sphere) {
            mini = primitive.col(0) - primitive.col(1);
            maxi = primitive.col(0) + primitive.col(1);
        } else {
            mini = primitive.rowwise().minCoeff();      // rowwise --- get coeff in each row
            maxi = primitive.rowwise().maxCoeff();
            Eigen::Vector3<Ty> diff = maxi - mini;
            for (int i = 0; i < 3; i++) {               // expand ill-posed AABB
                if (diff(i) < 1e-4) {
                    mini(i) -= 1e-4;
                    maxi(i) += 1e-4;
                }
            }
        }
    }

    Ty axis_centroid(SplitAxis axis) const {
        return (maxi[axis] + mini[axis]) * 0.5f;
    }

    AABB operator+(const AABB& aabb) const {
        return AABB(aabb.mini.cwiseMax(mini), aabb.maxi.cwiseMax(maxi));
    }

    AABB& operator+=(const AABB& aabb) {
        this->mini = aabb.mini.cwiseMin(this->mini);
        this->maxi = aabb.maxi.cwiseMax(this->maxi);
        return *this;
    }

    void clear() {
        mini.fill(1e4);
        maxi.fill(-1e4);
    }

    Ty area() const {
        Eigen::Vector3<Ty> diff = maxi - mini;
        return 2. * (diff(0) * diff(1) + diff(1) * diff(2) + diff(0) * diff(2));
    }
};

template <typename Ty>
struct BVHInfo {
    // BVH is for both triangle meshes and spheres
    AABB<Ty> bound;
    Eigen::Vector3<Ty> centroid;
    int prim_idx;
    int obj_idx;

    BVHInfo(): centroid(Eigen::Vector3<Ty>::Zero()), prim_idx(-1), obj_idx(-1) {}
    BVHInfo(const Eigen::Matrix3<Ty>& primitive, int prim_idx, int obj_idx, bool is_sphere = false): 
        bound(primitive, is_sphere), prim_idx(prim_idx), obj_idx(obj_idx) 
    {
        // Extract two vertices for the primitive, once converted to AABB
        // We don't need to distinguish between mesh or sphere
        // Note that vertex vectors in the primitive matrix are col-major order
        if (is_sphere)
            centroid = primitive.col(0);
        else
            centroid = primitive.rowwise().mean();      // barycenter
    }
};

template <typename Ty>
struct AxisBins {
    AABB<Ty> bound;
    int prim_cnt;

    AxisBins(): prim_cnt(0.f) {}

    void push(const BVHInfo<Ty>& bvh) {
        bound += bvh.bound;
        prim_cnt ++;
    }

    void push(const AABB<Ty>& aabb) {
        bound += aabb;
        prim_cnt ++;
    }
};

template <typename Ty>
class BVHNode {
public:
    BVHNode(): base(0), prim_num(0), axis(NONE), lchild(nullptr), rchild(nullptr) {}
    BVHNode(int base, int prim_num): base(base), prim_num(prim_num), axis(NONE), lchild(nullptr), rchild(nullptr) {}
    ~BVHNode() {
        if (lchild != nullptr) delete lchild;
        if (rchild != nullptr) delete rchild;
    }
public:
    // The axis start and end are scaled up a little bit
    SplitAxis max_extent_axis(const std::vector<BVHInfo<Ty>>& bvhs, std::vector<Ty>& bins) const;
public:
    int base;
    int prim_num;

    AABB<Ty> bound;
    SplitAxis axis;
    BVHNode<Ty> *lchild, *rchild;
};

/**
 * @note @attention
 * For taichi lang, there is actually a huge problem: 
 * To traverse the BVH tree, it is best to choose between lchild and rchild 
 * according to split axis and the current ray direction, since for example
 * if ray points along neg-x-axis, and the node is split along x axis
 * lchild contains smaller x-coordinate nodes and the otherwise for rchild
 * we should of course first traverse the rchild. BUT, taichi lang has neither
 * (1) kernel stack nor (2) thread local dynamic memory allocation, therefore
 * neither can we record whether the node is accessed nor store the node to be accessed
 * in a stack. Allocate a global field which might incur (H * W * D * 4) B memory consumption,
 * where D is the depth of BVH tree, which can be extremely inbalanced when the primitives
 * are distributed with poor uniformity. Therefore, in this implementation, we can only opt
 * for a suboptimal solution, to traverse the tree using just DFS, resulting in a simple
 * traversal method: the index for the linear node container is monotonously increasing. So 
 * in this implementation, split axis will not be included (which will not be used anyway).
 * TODO: we can opt for dynamic snode, but I think this would be ugly. Simple unidirectional traversal
 * would be much better than brute-force traversal with only AABB pruning.
 */
template <typename Ty>
class LinearNode {
public:
    LinearNode(): base(0), prim_cnt(0) {
        mini.resize({3});
        maxi.resize({3});
    }
    LinearNode(const BVHNode<Ty> *const bvh): base(bvh->base), prim_cnt(bvh->prim_num) {
        mini.resize({3});
        maxi.resize({3});
        Ty *const min_ptr = mini.mutable_data(0), *const max_ptr = maxi.mutable_data(0);
        for (int i = 0; i < 3; i++) {
            min_ptr[i] = bvh->bound.mini(i);
            max_ptr[i] = bvh->bound.maxi(i);
        }
    };       // linear nodes are initialized during DFS binary tree traversal
public:
    // The linearized BVH tree should contain: bound, base, prim_cnt, rchild_offset, total_offset (to skip the entire node)
    py::array_t<Ty> mini;
    py::array_t<Ty> maxi;
    int base, prim_cnt;         // indicate the starting point and the length of the node
    int all_offset;             // indicate the rchild pos and the offset to the next node
};

template <typename Ty>
class LinearBVH {
public:
    LinearBVH(): obj_idx(-1), prim_idx(-1) {
        mini.resize({3});
        maxi.resize({3});
    }
    LinearBVH(const BVHInfo<Ty>& bvh): obj_idx(bvh.obj_idx), prim_idx(bvh.prim_idx) {
        mini.resize({3});
        maxi.resize({3});
        Ty *const min_ptr = mini.mutable_data(0), *const max_ptr = maxi.mutable_data(0);
        for (int i = 0; i < 3; i++) {
            min_ptr[i] = bvh.bound.mini(i);
            max_ptr[i] = bvh.bound.maxi(i);
        }
    }
public:
    int obj_idx;
    int prim_idx;
    py::array_t<Ty> mini;
    py::array_t<Ty> maxi;
};

}; // end name space scds
