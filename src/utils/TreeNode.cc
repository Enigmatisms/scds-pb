#include "utils/TreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim, size_t Nchild>
std::shared_ptr<TreeNode<Ty, Ndim, Nchild>> TreeNode<Ty, Ndim, Nchild>::try_get_child(size_t child_idx) {
    ProfilePhase _(Prof::TreeNodeTryGetChild);
    if (!childs[child_idx]) {
        // if possible, we add another depth
        Pointx half_size = size * 0.5;
        Pointx offset = get_child_offset(half_size, child_idx);
        {
            ProfilePhase _(Prof::TreeNodeCreateSharedPtr);
            
            childs[child_idx] = std::make_shared<TreeNode>(center + offset, half_size, this->shared_from_this());
        }
        #ifdef TREE_NODE_MEMORY_PROFILE
            sharedPtrBytes += sizeof(std::shared_ptr<TreeNode>);
        #endif //TREE_NODE_MEMORY_PROFILE
    }
    return childs[child_idx];
}

template<typename Ty, size_t Ndim, size_t Nchild>
Point<Ty, Ndim> TreeNode<Ty, Ndim, Nchild>::get_child_offset(const Pointx& half_size, size_t child_id) {
    ProfilePhase _(Prof::TreeNodeGetChildOffset);

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

template<typename Ty, size_t Ndim, size_t Nchild>
size_t TreeNode<Ty, Ndim, Nchild>::get_size() const {
    size_t size = sizeof(*this);  // Size of the current instance
    size += sub_idxs.size() * sizeof(int);
    for (const auto& child : childs) {
        if (child)
            size += child->get_size();
    }
    return size;
}

template class TreeNode<float, 2, 4>;               // 2D Quad-tree (vector container)
template class TreeNode<float, 3, 8>;               // 3D Octree    (vector container)
template class TreeNode<double, 2, 4>;              // 2D Quad-tree (vector container)
template class TreeNode<double, 3, 8>;              // 3D Octree    (vector container) 

} // namespace name
