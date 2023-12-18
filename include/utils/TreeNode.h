#include <array>
#include <vector>
#include <memory>
#include <unordered_set>
#include "Point.h"
namespace scds {

template<typename T, size_t Ndim, size_t Nchild>
class TreeNode {
using Pointx = Point<T, Ndim>;
public:
    TreeNode(
        const Pointx& tl, const Pointx& br, 
        std::unordered_set<size_t>&& idxs,
        std::shared_ptr<TreeNode> parent, 
        std::shared_ptr<std::vector<Pointx>> pts_ptr
    );

    TreeNode(
        const Pointx& center, const Pointx& size, 
        std::unordered_set<size_t>&& idxs,
        std::shared_ptr<TreeNode> parent, 
        std::shared_ptr<std::vector<Pointx>> pts_ptr
    );
public:
    bool point_in_range(const Pointx& pt) const {
        auto diff = (center - pt).abs();
        return ((size - diff) > 0).all();
    }

    bool exist(size_t idx) const {
        return sub_idxs.count(idx);
    }

    bool remove(size_t idx) const {
        auto it = sub_idxs.find(idx);
        if (it != sub_idxs.end()) {
            sub_idxs.erase(it);
            return true;
        }
        return false;                   // not removed since the target does not exist
    }

    auto get_parent() const {
        if (auto ptr = parent.lock())
            return ptr;
        else
            return nullptr;
    }
protected:
    // pointer to the parent node
    std::weak_ptr<TreeNode> parent;

    // pointer to child nodes
    std::array<std::shared_ptr<TreeNode>, Nchild> childs;

    // center to the sub-tree (partitioned space)
    Pointx center;

    // extent of the partitioned space
    Pointx size;

    // whether the node is a leaf node
    bool is_leaf;

    // global chunk of memory that is visible to all tree nodes.
    std::shared_ptr<std::vector<Pointx>> pts;

    // the indices to points contained in this sub-tree
    std::unordered_set<size_t> sub_idxs;
    // question is: where should we store the nodes and how to recursive free the tree 
};

// node for KD-tree
template<typename T, size_t Ndim>
using KdNode = TreeNode<T, Ndim, 2>;        
template<typename T>
// node for Quad-Tree
using QdNode = TreeNode<T, 2, 4>;
// node for Octree
template<typename T>
using OcNode = TreeNode<T, 3, 8>;

} // end namespace scds