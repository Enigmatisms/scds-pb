#pragma once
/**
 * KD-tree for primitives. This will be migrated to my AdaPT renderer
 * @author: Qianyue He
 * @date:   2024-2-22
*/

#include <array>
#include <vector>
#include <memory>
#include <Eigen/Core>
#include "Point.h"
#include "utils/stats.h"
#include "bvh/bvh_helper.h"

namespace scds {

/**
 * @brief KD-tree efficient nearest neighbor search structure
*/

#define KD_PRIM_NODE_MEMORY_PROFILE
#ifdef KD_PRIM_NODE_MEMORY_PROFILE
STAT_MEMORY_COUNTER("PrimNode/Total TreeNode Size", primNodeBytes);
STAT_COUNTER("PrimNode/Node count", primNodeCount);
STAT_COUNTER("PrimNode/Leaf nodes", primLeafNodes);
#endif // KD_PRIM_NODE_MEMORY_PROFILE

template<typename Ty>
class KDPrimNode {
using This = KDPrimNode<Ty>;
using Vec3 = Eigen::Vector3<Ty>;
using Mat3 = Eigen::Matrix3<Ty>;
public:
    template <typename AABBType>
    KDPrimNode(
        AABBType&& aabb, std::vector<int>&& idxs,
        SplitAxis axis = SplitAxis::NONE, Ty split_pos = 0 
    ): aabb(std::forward<AABBType>(aabb)),
        split_axis(axis), split_pos(split_pos), 
        sub_idxs(std::make_unique<std::vector<int>>(std::move(idxs)))
    {
        num_elems = static_cast<int>(sub_idxs->size());
        #ifdef KD_PRIM_NODE_MEMORY_PROFILE
            primNodeBytes += this->get_size();
            primLeafNodes ++;
            primNodeCount ++;
        #endif //KD_PRIM_NODE_MEMORY_PROFILE
    }

    template <typename Ptype1, typename Ptype2>
    KDPrimNode(
        Ptype1&& mini, Ptype2&& maxi, 
        std::vector<int>&& idxs,
        SplitAxis axis = SplitAxis::NONE, Ty split_pos = 0 
    ): aabb(std::forward<Ptype1>(mini), std::forward<Ptype2>(maxi)), 
        split_axis(axis), split_pos(split_pos),
        sub_idxs(std::make_unique<std::vector<int>>(std::move(idxs)))
    {
        num_elems = static_cast<int>(sub_idxs->size());
        #ifdef KD_PRIM_NODE_MEMORY_PROFILE
            primNodeCount ++;
            primNodeBytes += this->get_size();
            primLeafNodes ++;
        #endif //KD_PRIM_NODE_MEMORY_PROFILE
    }
public:
    void insert(int value) {
        sub_idxs->emplace_back(value);
        #ifdef KD_PRIM_NODE_MEMORY_PROFILE
            primNodeBytes += sizeof(int);
        #endif //KD_PRIM_NODE_MEMORY_PROFILE
    }
    // query child node with given index. If the queried child is nullptr, return directly.

    void split_leaf_node(
        const std::vector<AABB<Ty>>& aabbs, Ty split_pos, SplitAxis split_axis,
        AABB<Ty>&& lchild_aabb, AABB<Ty>&& rchild_aabb
    );
    
    std::shared_ptr<KDPrimNode> lchild() {
        return _lchild;
    }

    std::shared_ptr<KDPrimNode> rchild() {
        return _rchild;
    }

    // get the indices stored in the node
    const std::vector<int>& get_indices() const {
        return *sub_idxs;
    }

    // get the orthogonal axis which span the maximum extent
    SplitAxis max_extent_axis(const std::vector<AABB<Ty>>& aabbs, std::vector<Ty>& bins, std::vector<AxisBins<Ty>>& idx_bins) const;

    template <typename AABBType>
    auto add_child(
        AABBType&& aabb, std::vector<int>&& pt_idxs, 
        bool is_left = false, SplitAxis _split_axis = SplitAxis::NONE, Ty _split_pos = 0 
    ) {
        ProfilePhase _(Prof::TreeNodeAddChild);
        {
            auto new_child = std::make_shared<KDPrimNode>(
                std::forward<Ptype1>(aabb), 
                std::move(pt_idxs), 
                _split_axis, _split_pos
            );
            if (is_left)
                _lchild = std::move(new_child);
            else
                _rchild = std::move(new_child);
        }
    }

    /**
     * @brief Get the child node that is closer to the current query point
     * @param query:       top right corner of the search range (size + radius)
    */ 
    template <typename PointType>
    std::shared_ptr<KDPrimNode> nearest_child(PointType&& query) const {
        Ty diff_l = (_lchild->center - query).length2(),
           diff_r = (_rchild->center - query).length2();
        return diff_l <= diff_r ? _lchild : _rchild;
    }

    template <typename PointType>
    bool in_range(PointType&& query) const {
        for (int i = 0; i < 3; i++) {
            if (query[i] < aabb.min[i]) return false;
            else if (query[i] > aabb.max[i]) return false;
        }
        return true;
    }

    template <typename PointType>
    std::shared_ptr<KDPrimNode> resident_child(PointType&& query) const {
        if (_lchild->in_range(std::forward<PointType>(query)))
            return _lchild;
        else if (_rchild->in_range(std::forward<PointType>(query)))
            return _rchild;
        return nearest_child(std::forward<PointType>(query));
    }

    std::shared_ptr<KDPrimNode> the_other(const KDPrimNode* const node_ptr) const {
        if (_lchild.get() == node_ptr)
            return _rchild;
        else if (_rchild.get() == node_ptr)
            return _lchild;
        std::cerr << "node_ptr is not the child ptr of the current node." << std::endl;
        return {nullptr};
    }

    bool is_leaf() const noexcept {
        return split_axis == SplitAxis::NONE;
    }

    void set_leaf() noexcept {
        // make the current node a leaf node (when sub_idxs is not a nullptr)
        sub_idxs = std::make_unique<std::vector<int>>();
        split_axis = SplitAxis::NONE;
        _lchild.reset();
        _rchild.reset();
    }

    void set_non_leaf() noexcept {
        // make the current node a non-leaf node (when sub_idxs is a nullptr)
        sub_idxs.reset(nullptr);
    }
public:
    AABB<Ty> aabb;

    // Splitting axis
    SplitAxis split_axis;
    Ty split_pos;

    // number of points stored in the subtree
    int num_elems;
    std::unique_ptr<std::vector<int>> sub_idxs;
private:
    size_t get_size() const;
protected:
    std::shared_ptr<KDPrimNode> _lchild;
    std::shared_ptr<KDPrimNode> _rchild;
    // indices will only be stored in the leaf nodes (for non-leaf nodes, this is a nullptr)
};

template <typename Ty>
struct LinearKdNode {
    LinearKdNode(): idx_base(0), idx_base(0), r_offset(0), all_offset(0) {
        mini.resize({3});
        maxi.resize({3});
    }
    LinearKdNode(const KDPrimNode<Ty> *const kd_node): 
        idx_base(0), idx_base(0), r_offset(0), all_offset(0) 
    {
        mini.resize({3});
        maxi.resize({3});
        Ty *const min_ptr = mini.mutable_data(0), *const max_ptr = maxi.mutable_data(0);
        for (int i = 0; i < 3; i++) {
            min_ptr[i] = kd_node->aabb.mini(i);
            max_ptr[i] = kd_node->aabb.maxi(i);
        }
    };       // linear nodes are initialized during DFS binary tree traversal

    py::array_t<float> mini;
    py::array_t<float> maxi;
    int idx_base, idx_num;          // starting index of the linearized sub_idxs, number of indices
    int r_offset, all_offset;       // offset to the rchild and the number of nodes in this (sub-)tree
};

/**
 * 
*/
/**
 * we only return the root node (linearize and export to Taichi end)
 * main logic walk-through
 * 
 * special note:
 * (1) the tree construction will be a recursive process: primitives are indexed with indices
 * for the lchild and rchild of a node, should we re-order the primitives to accelerate querying (cache)?
 * I don't think so, since one primitive can exist in different nodes, this re-ordering might 
 * be tricky? Maybe it is only feasible for BVH.
 * (2) the splitting plane is chosen to be one of the planes of the bounding box
 * (3) the tree nodes will be linearized in a pre-order way
 * (4) the primitives fall within a certain leave node can be extremely unbalanced, since we are splitting the tree
 * according to SAH: small primitive clusters will result in leaf nodes with more primitives and vice-versa
 * So we might need **dynamic SNode in Taichi lang**
 * 
 * special note for ray intersection logic:
 * (1) Note that, before this is used in Taichi end, node linearization will be performed (pre-order)
 * We should be clear whether pre-order linearization or in-order linearization is better (since PBR-book performs
 * in-order linearization). Here we assume the linearization is pre-order one, therefore (root)->(l-tree)->(r-tree).
 * 
 * (2) This runs into the same pit-fall: in Taichi it will be un-easy to do optimal ray-intersection test, since
 * stack (todo-list) will not be easy to implement, unless we are using tile-based rendering (no-preview): pre-allocated
 * stack (with fixed upper limit of stack capacity) can be reused between different tiles. Otherwise, we should find a way
 * to allocate a local / variable array. Maybe I don't have to? I can simply allocate a fixed size stack (vec16 will be ok)
 * depending on the user setting, we can adjust the vec length (for example, 'coarse' --- vec4, 'medium' --- vec8, 'fine' --- vec16)
 * 
 * ray intersection logic:
 * 1. intersect the current node, check whether the ray intersects the AABB or not, if not: skip the whole node (we store an offset
 * in the node)
 * 2. for pre-order organization, we store a offset to the rchild, so here we first consider which child to examine first, depending on
 * whether the ray intersects the lchild or rchild first. The other one will be cached inside the local stack (along with AABB t_min)
 * if the traverse of one leaf-node yields a intersection t smaller than the t_min recorded in the stack, the top of the stack can
 * be directly popped
 * 
 * Objective: implement two versions of ray-intersection (one in Taichi, one in C++, due 2.25)
 * Other parts are yet to be cleared
*/
}       // end namespace scds