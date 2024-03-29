#pragma once
#include <array>
#include <vector>
#include <memory>
#include "Point.h"
#include "utils/stats.h"
#include "utils/utils.h"

namespace scds {

/**
 * k-d tree for points
 * A better tree representation for binary tree is possible (linearized)
*/

#define BIN_TREE_NODE_MEMORY_PROFILE
#ifdef BIN_TREE_NODE_MEMORY_PROFILE
STAT_MEMORY_COUNTER("BinTreeNode/Total TreeNode Size", binTreeBytes);
STAT_COUNTER("BinTreeNode/Node count", binNodeCount);
STAT_COUNTER("BinTreeNode/Leaf nodes", binLeafNodes);
#endif // TREE_NOD_MEMORY_PROFILE

template<typename Ty, size_t Ndim>
class BinTreeNode {
using Pointx = Point<Ty, Ndim>;
using This   = BinTreeNode<Ty, Ndim>;
public:
    /**
     * Note: you may wonder why I use two template types here to have two universal references
     * I think (this is my current understanding - 12.20), if center and size are of different types (l/rvalue)
     * one unified template type might lead to compilation error
    */
    template <typename Ptype1, typename Ptype2>
    BinTreeNode(
        Ptype1&& center, Ptype2&& size, 
        SplitAxis axis = SplitAxis::NONE, Ty split_pos = 0 
    ): center(std::forward<Ptype1>(center)), half_size(std::forward<Ptype2>(size)),
        split_axis(axis), split_pos(split_pos), num_elems(0) 
    {
        sub_idxs = std::make_unique<std::vector<int>>();
        #ifdef BIN_TREE_NODE_MEMORY_PROFILE
            binTreeBytes += this->get_size();
            binLeafNodes ++;
            binNodeCount ++;
        #endif //BIN_TREE_NODE_MEMORY_PROFILE
    }

    template <typename Ptype1, typename Ptype2>
    BinTreeNode(
        Ptype1&& center, Ptype2&& size, 
        std::vector<int>&& idxs,
        SplitAxis axis = SplitAxis::NONE, Ty split_pos = 0 
    ): center(std::forward<Ptype1>(center)), half_size(std::forward<Ptype2>(size)), 
        split_axis(axis), split_pos(split_pos),
        sub_idxs(std::make_unique<std::vector<int>>(std::move(idxs)))
    {
        num_elems = static_cast<int>(sub_idxs->size());
        #ifdef BIN_TREE_NODE_MEMORY_PROFILE
            binNodeCount ++;
            binTreeBytes += this->get_size();
            binLeafNodes ++;
        #endif //BIN_TREE_NODE_MEMORY_PROFILE
    }
public:
    void insert(int value) {
        sub_idxs->emplace_back(value);
        #ifdef BIN_TREE_NODE_MEMORY_PROFILE
            binTreeBytes += sizeof(int);
        #endif //BIN_TREE_NODE_MEMORY_PROFILE
    }
    // query child node with given index. If the queried child is nullptr, return directly.

    void split_leaf_node(const std::vector<Pointx>& pts);
    
    std::shared_ptr<BinTreeNode> lchild() {
        return _lchild;
    }

    std::shared_ptr<BinTreeNode> rchild() {
        return _rchild;
    }

    // get the indices stored in the node
    const std::vector<int>& get_indices() const {
        return *sub_idxs;
    }

    // get the orthogonal axis which span the maximum extent
    SplitAxis max_extent_axis(const std::vector<Pointx>& pts) const;

    template <typename Ptype1, typename Ptype2>
    auto add_child(
        Ptype1&& ctr, Ptype2&& new_size, std::vector<int>&& pt_idxs, 
        bool is_left = false, SplitAxis _split_axis = SplitAxis::NONE, Ty _split_pos = 0 
    ) {
        ProfilePhase _(Prof::TreeNodeAddChild);
        {
            auto new_child = std::make_shared<BinTreeNode>(
                std::forward<Ptype1>(ctr), 
                std::forward<Ptype2>(new_size), 
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
    std::shared_ptr<BinTreeNode> nearest_child(PointType&& query) const {
        Ty diff_l = (_lchild->center - query).length2(),
           diff_r = (_rchild->center - query).length2();
        return diff_l <= diff_r ? _lchild : _rchild;
    }

    template <typename PointType>
    bool in_range(PointType&& query) const {
        return (query >= center - half_size).all() && (query < center + half_size).all();
    }

    template <typename PointType>
    std::shared_ptr<BinTreeNode> resident_child(PointType&& query) const {
        if (_lchild->in_range(std::forward<PointType>(query)))
            return _lchild;
        else if (_rchild->in_range(std::forward<PointType>(query)))
            return _rchild;
        return nearest_child(std::forward<PointType>(query));
    }

    std::shared_ptr<BinTreeNode> the_other(const BinTreeNode* const node_ptr) const {
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
    Pointx center;
    Pointx half_size;

    // Splitting axis
    SplitAxis split_axis;
    Ty split_pos;

    // number of points stored in the subtree
    int num_elems;
    std::unique_ptr<std::vector<int>> sub_idxs;
private:
    size_t get_size() const;
protected:
    std::shared_ptr<BinTreeNode> _lchild;
    std::shared_ptr<BinTreeNode> _rchild;
    // indices will only be stored in the leaf nodes (for non-leaf nodes, this is a nullptr)
};

template<typename T>
// node for Quad-Tree
using KdNode2 = BinTreeNode<T, 2>;
// node for Octree
template<typename T>
using KdNode  = BinTreeNode<T, 3>;

} // end namespace scds