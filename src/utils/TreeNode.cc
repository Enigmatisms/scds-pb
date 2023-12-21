#include "utils/TreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim, size_t Nchild>
std::shared_ptr<TreeNode<Ty, Ndim, Nchild>> TreeNode<Ty, Ndim, Nchild>::try_get_child(size_t child_idx) {
    if (!childs[child_idx]) {
        // if possible, we add another depth
        Pointx half_size = size / 2;
        Pointx offset = get_child_offset(half_size, child_idx);
        childs[child_idx] = std::make_shared<TreeNode>(center + offset, half_size, this->shared_from_this(), pts);
    }
    return childs[child_idx];
}

template class TreeNode<float, 2, 4>;               // 2D Quad-tree
template class TreeNode<float, 3, 8>;               // 3D Octree
template class TreeNode<double, 2, 4>;              // 2D Quad-tree
template class TreeNode<double, 3, 8>;              // 3D Octree

} // namespace name
