#include "utils/TreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim, size_t Nchild>
std::shared_ptr<TreeNode<Ty, Ndim, Nchild>> TreeNode<Ty, Ndim, Nchild>::try_get_child(size_t child_idx) {
    ProfilePhase _(Prof::TreeNodeTryGetChild);
    if (!childs[child_idx]) {
        // if possible, we add another depth
        Pointx half_size = size / 2;
        Pointx offset = get_child_offset(half_size, child_idx);
        childs[child_idx] = std::make_shared<TreeNode>(center + offset, half_size, this->shared_from_this(), pts);
    }
    return childs[child_idx];
}

template<typename Ty, size_t Ndim, size_t Nchild>
size_t TreeNode<Ty, Ndim, Nchild>::get_size() const {
    size_t size = sizeof(*this);  // Size of the current instance

    // Add the size of the unordered_set
    size += sub_idxs.bucket_count() * sizeof(size_t);

    // Add the sizes of child nodes
    for (const auto& child : childs) {
        if (child)
            size += child->get_size();
    }

    return size;
}

template class TreeNode<float, 2, 4>;               // 2D Quad-tree
template class TreeNode<float, 3, 8>;               // 3D Octree
template class TreeNode<double, 2, 4>;              // 2D Quad-tree
template class TreeNode<double, 3, 8>;              // 3D Octree

} // namespace name
