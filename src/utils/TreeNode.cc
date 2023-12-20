#include "utils/TreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim, size_t Nchild>
template <typename Ptype1, typename Ptype2>
TreeNode<Ty, Ndim, Nchild>::TreeNode(
    Ptype1&& center, Ptype2&& size, 
    std::weak_ptr<This> parent, 
    std::shared_ptr<std::vector<Pointx>> pts_ptr
): center(std::forward<Ptype1>(center)), size(std::forward<Ptype1>(size)), parent(parent), pts(pts_ptr)  {
    // TODO:
}

template<typename Ty, size_t Ndim, size_t Nchild>
template <typename Ptype1, typename Ptype2>
TreeNode<Ty, Ndim, Nchild>::TreeNode(
    Ptype1&& center, Ptype2&& size, 
    std::weak_ptr<This> parent, 
    std::unordered_set<size_t>&& sub_idxs,
    std::shared_ptr<std::vector<Pointx>> pts_ptr
):  center(std::forward<Ptype1>(center)), size(std::forward<Ptype1>(size)), 
    parent(parent), pts(pts_ptr), sub_idxs(std::move(sub_idxs))  
{
    // TODO: 
}

template<typename Ty, size_t Ndim, size_t Nchild>
auto TreeNode<Ty, Ndim, Nchild>::get_child(size_t child_idx) {
    if (childs[child_idx] != nullptr) {
        return child_idx;
    } else {
        // if possible, we add another depth
        Pointx half_size = size / 2;
        Pointx offset = get_child_offset(half_size, child_idx);
        childs[child_idx] = std::make_shared<TreeNode>(center + offset, half_size, shared_from_this(), pts);
        return childs[child_idx];
    }
}


template<typename Ty, size_t Ndim, size_t Nchild>
template <typename PointType>
Point<Ty, Ndim> TreeNode<Ty, Ndim, Nchild>::get_child_offset(PointType&& half_size, size_t child_id) {
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

template class TreeNode<float, 2, 4>;               // 2D Quad-tree
template class TreeNode<float, 3, 8>;               // 3D Octree
template class TreeNode<double, 2, 4>;              // 2D Quad-tree
template class TreeNode<double, 3, 8>;              // 3D Octree

} // namespace name
