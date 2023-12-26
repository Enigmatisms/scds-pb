#pragma once
#include <array>
#include <vector>
#include <memory>
#include <unordered_set>
#include "Point.h"
namespace scds {

/**
 * Notes from PBR-book (4.4 Kd-Tree): about (1) save memory (2) continous memory block can improve cache performance
 * > Rather than storing two pointers or offsets, we lay the nodes out in a way that lets us only store one child pointer:
 * all of the nodes are allocated in a single contiguous block of memory, and the child of an interior node that is
 * responsible for space below the splitting plane is always stored in the array position immediately after its parent
 * (this layout also improves cache performance, by keeping at least one child close to its parent in memory)
 * 
 * I did this before, for BVH tree (linearize the tree, since they are binary, easy to derive)
 * I suspect that here we can do the same optimization, but I will stick to the simpler version first
 * Then I might go back and try to improve this
*/

template<typename Ty, size_t Ndim, size_t Nchild>
class TreeNode : public std::enable_shared_from_this<TreeNode<Ty, Ndim, Nchild>> {
using Pointx = Point<Ty, Ndim>;
using This   = TreeNode<Ty, Ndim, Nchild>;
public:
    /**
     * Note: you may wonder why I use two template types here to have two universal references
     * I think (this is my current understanding - 12.20), if center and size are of different types (l/rvalue)
     * one unified template type might lead to compilation error
    */
    template <typename Ptype1, typename Ptype2>
    TreeNode(
        Ptype1&& center, Ptype2&& size, 
        std::weak_ptr<TreeNode> parent, 
        std::shared_ptr<std::vector<Pointx>> pts_ptr
    ): center(std::forward<Ptype1>(center)), size(std::forward<Ptype2>(size)), parent(parent), pts(pts_ptr)  {}

    template <typename Ptype1, typename Ptype2>
    TreeNode(
        Ptype1&& center, Ptype2&& size, 
        std::weak_ptr<TreeNode> parent, 
        std::unordered_set<size_t>&& sub_idxs,
        std::shared_ptr<std::vector<Pointx>> pts_ptr
    ):  center(std::forward<Ptype1>(center)), size(std::forward<Ptype2>(size)), 
    parent(parent), pts(pts_ptr), sub_idxs(std::move(sub_idxs)) {}

    template <typename PointType>
    static Pointx get_child_offset(PointType&& half_size, size_t child_id) {
        Pointx offset;
        for (size_t i = 0; i < Ndim; i++) {
            size_t index = Ndim - 1 - i;
            if (child_id & 1)
                offset[index] = half_size[index];
            else
                offset[index] = -half_size[index];
            child_id >>= 1;
        }
        return offset;
    }
public:
    // query child node with given index. If the queried child is nullptr, we will create a new child node inplace and return it
    std::shared_ptr<TreeNode<Ty, Ndim, Nchild>> try_get_child(size_t child_idx);

    void insert(size_t index) {
        sub_idxs.emplace(index);
    }

    // query child node with given index. If the queried child is nullptr, return directly.
    std::shared_ptr<TreeNode> get_child(size_t child_idx) {
        return childs[child_idx];
    }

    // get the indices stored in the node
    const std::unordered_set<size_t>& get_indices() const {
        return sub_idxs;
    }

    bool point_in_range(const Pointx& pt) const {
        auto diff = (center - pt).abs();
        return ((size - diff) > 0).all();
    }

    bool exist(size_t idx) const {
        return sub_idxs.count(idx);
    }

    bool remove(size_t idx) {
        auto it = sub_idxs.find(idx);
        if (it != sub_idxs.end()) {
            sub_idxs.erase(it);
            return true;
        }
        return false;                   // not removed since the target does not exist
    }

    std::shared_ptr<TreeNode> get_parent() const {
        if (auto ptr = parent.lock())
            return ptr;
        else
            return {};
    }

    size_t num_points() const {
        return sub_idxs.size();
    }

    template <typename Ptype1, typename Ptype2>
    auto add_child(Ptype1&& ctr, Ptype2&& new_size, std::unordered_set<size_t>&& pt_idxs, size_t child_id) {
        childs[child_id] = std::make_shared<TreeNode>(
            std::forward<Ptype1>(ctr), 
            std::forward<Ptype2>(new_size), 
            this->shared_from_this(), std::move(pt_idxs), this->pts
        );
    }

    /**
     * @brief Check if a node overlaps with a specified search range
     * @param tr:       top right corner of the search range (size + radius)
     * @param bl:       bottom left corner of the search range (size - radius)
     * @param node_ptr: the node to check
    */ 
    template <typename PointType>
    inline bool overlap_range(PointType&& tr, PointType&& bl) const {
        auto tr_node = this->center + this->size;
        auto bl_node = this->center - this->size;
        for (size_t i = 0; i < Ndim; i++) {
            if (tr_node[i] < bl[i] || tr[i] < bl_node[i]) return false;
        }
        return true;
    }
public:
    // center to the sub-tree (partitioned space)
    Pointx center;

    // extent of the partitioned space
    Pointx size;

    // whether the node is a leaf node
    bool is_leaf;
protected:
    // pointer to the parent node
    std::weak_ptr<TreeNode> parent;

    // pointer to child nodes
    std::array<std::shared_ptr<TreeNode>, Nchild> childs;

    // global chunk of memory that is visible to all tree nodes.
    std::shared_ptr<std::vector<Pointx>> pts;

    // the indices to points contained in this sub-tree
    std::unordered_set<size_t> sub_idxs;
    // question is: where should we store the nodes and how to recursive free the tree 
};

template<typename T>
// node for Quad-Tree
using QdNode = TreeNode<T, 2, 4>;
// node for Octree
template<typename T>
using OcNode = TreeNode<T, 3, 8>; 

} // end namespace scds