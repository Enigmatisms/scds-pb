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
class TreeNode : public std::enable_shared_from_this<TreeNode> {
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
    );

    TreeNode(
        Ptype1&& center, Ptype2&& size, 
        std::weak_ptr<TreeNode> parent, 
        std::unordered_set<size_t>&& sub_idxs,
        std::shared_ptr<std::vector<Pointx>> pts_ptr
    );
public:
    auto get_child(size_t child_idx);

    /// You need to consider both the following situation: (1) the tree structure is not yet built (2) insert in an existing node
    // Add a batch of point indices (to form a child tree)
    template <typename Ptype1, typename Ptype2>
    auto add_child(Ptype1&& ctr, Ptype2&& new_size, std::unordered_set<size_t>&& pt_idxs, size_t child_id) {
        childs[child_id] = std::make_shared<TreeNode>(
            std::forward<Ptype1>(ctr), 
            std::forward<Ptype2>(new_size), 
            shared_from_this(), std::move(pt_idxs), this->pts
        );
    }

    template <typename PointType>
    static Pointx quadrant_offset(PointType&& half_size, size_t quad_id);

    // Add one point index (to form a child tree)

    void insert(size_t index) {
        sub_idxs.emplace(index);
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