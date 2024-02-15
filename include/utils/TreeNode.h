#pragma once
#include <array>
#include <vector>
#include <memory>
#include "Point.h"
#include "utils/stats.h"

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

#define TREE_NODE_MEMORY_PROFILE
#ifdef TREE_NODE_MEMORY_PROFILE
STAT_MEMORY_COUNTER("TreeNode/Total TreeNode Size", treeBytes);
STAT_MEMORY_COUNTER("TreeNode/Shared Ptr Size", sharedPtrBytes);
STAT_COUNTER("TreeNode/Node count", nodeCount);
STAT_COUNTER("TreeNode/Leaf nodes", leafNodes);
#endif // TREE_NOD_MEMORY_PROFILE

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
        std::weak_ptr<TreeNode> parent
    ): center(std::forward<Ptype1>(center)), size(std::forward<Ptype2>(size)), num_points(0), parent(parent) {
        sub_idxs = std::make_unique<std::vector<int>>();
        #ifdef TREE_NODE_MEMORY_PROFILE
            treeBytes += this->get_size();
            sharedPtrBytes += sizeof(std::shared_ptr<std::vector<Pointx>>);
            leafNodes ++;
            nodeCount ++;
        #endif //TREE_NODE_MEMORY_PROFILE
    }

    template <typename Ptype1, typename Ptype2>
    TreeNode(
        Ptype1&& center, Ptype2&& size, 
        std::weak_ptr<TreeNode> parent, 
        std::vector<int>&& idxs
    ): center(std::forward<Ptype1>(center)), size(std::forward<Ptype2>(size)), 
    parent(parent), sub_idxs(std::make_unique<std::vector<int>>(std::move(idxs))) {
        num_points = static_cast<int>(sub_idxs->size());
        #ifdef TREE_NODE_MEMORY_PROFILE
            nodeCount ++;
            treeBytes += this->get_size();
            sharedPtrBytes += sizeof(std::shared_ptr<std::vector<Pointx>>);
            leafNodes ++;
        #endif //TREE_NODE_MEMORY_PROFILE
    }

    static Pointx get_child_offset(const Pointx& half_size, size_t child_id);
public:
    // query child node with given index. If the queried child is nullptr, we will create a new child node inplace and return it
    std::shared_ptr<TreeNode<Ty, Ndim, Nchild>> try_get_child(size_t child_idx);

    void insert(int value) {
        sub_idxs->emplace_back(value);
        #ifdef TREE_NODE_MEMORY_PROFILE
            treeBytes += sizeof(int);
        #endif //TREE_NODE_MEMORY_PROFILE
    }
    // query child node with given index. If the queried child is nullptr, return directly.
    std::shared_ptr<TreeNode> get_child(size_t child_idx) {
        return childs[child_idx];
    }

    // get the indices stored in the node
    const std::vector<int>& get_indices() const {
        return *sub_idxs;
    }

    std::shared_ptr<TreeNode> get_parent() const {
        if (auto ptr = parent.lock())
            return ptr;
        else
            return {};
    }

    template <typename Ptype1, typename Ptype2>
    auto add_child(Ptype1&& ctr, Ptype2&& new_size, std::vector<int>&& pt_idxs, size_t child_id) {
        ProfilePhase _(Prof::TreeNodeAddChild);
        {
            ProfilePhase _(Prof::TreeNodeCreateSharedPtr);
            childs[child_id] = std::make_shared<TreeNode>(
                std::forward<Ptype1>(ctr), 
                std::forward<Ptype2>(new_size), 
                this->shared_from_this(), std::move(pt_idxs)
            );
        }
        #ifdef TREE_NODE_MEMORY_PROFILE
            sharedPtrBytes += sizeof(std::shared_ptr<TreeNode>);
        #endif //TREE_NODE_MEMORY_PROFILE
    }

    /**
     * @brief Check if a node overlaps with a specified search range
     * @param tr:       top right corner of the search range (size + radius)
     * @param bl:       bottom left corner of the search range (size - radius)
     * @param node_ptr: the node to check
    */ 
    template <typename PointType>
    inline bool overlap_range(PointType&& tr, PointType&& bl) const {
        ProfilePhase _(Prof::NumProfCategories);
        auto tr_node = this->center + this->size;
        auto bl_node = this->center - this->size;
        for (size_t i = 0; i < Ndim; i++) {
            if (tr_node[i] < bl[i] || tr[i] < bl_node[i]) return false;
        }
        return true;
    }

    void overwrite_sub_idxs(std::vector<int>&& src) {
        #ifdef TREE_NODE_MEMORY_PROFILE
            if (sub_idxs) {
                auto size_of_set = sizeof(this->sub_idxs->size() * sizeof(int)) + sizeof(sub_idxs);
                treeBytes -= size_of_set;
            }
            treeBytes += sizeof(src.size() * sizeof(int)) + sizeof(src);
        #endif //TREE_NODE_MEMORY_PROFILE
        sub_idxs = std::make_unique<std::vector<int>>(std::move(src));
        num_points = static_cast<int>(sub_idxs->size());
    }

    bool is_leaf() const noexcept {
        return (sub_idxs != nullptr);
    }

    void set_leaf() noexcept {
        // make the current node a leaf node (when sub_idxs is not a nullptr)
        sub_idxs = std::make_unique<std::vector<int>>();
    }

    void set_non_leaf() noexcept {
        // make the current node a non-leaf node (when sub_idxs is a nullptr)
        sub_idxs.reset(nullptr);
    }
public:
    // center to the sub-tree (partitioned space)
    Pointx center;

    // extent of the partitioned space
    Pointx size;

    // number of points stored in the subtree
    int num_points;
private:
    size_t get_size() const;
protected:
    // pointer to the parent node
    std::weak_ptr<TreeNode> parent;

    // pointer to child nodes
    std::array<std::shared_ptr<TreeNode>, Nchild> childs;

    // indices will only be stored in the leaf nodes (for non-leaf nodes, this is a nullptr)
    std::unique_ptr<std::vector<int>> sub_idxs;
};

template<typename T>
// node for Quad-Tree
using QdNode = TreeNode<T, 2, 4>;
// node for Octree
template<typename T>
using OcNode = TreeNode<T, 3, 8>;

} // end namespace scds