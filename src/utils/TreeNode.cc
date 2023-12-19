#include "utils/TreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim, size_t Nchild>
TreeNode<Ty, Ndim, Nchild>::TreeNode(
    const Pointx& tl, const Pointx& br, 
    std::unordered_set<size_t>&& idxs,
    std::shared_ptr<This> parent, 
    std::shared_ptr<std::vector<Pointx>> pts_ptr
) {
    
}

template<typename Ty, size_t Ndim, size_t Nchild>
template<typename T>
auto TreeNode<Ty, Ndim, Nchild>::insert(T&& pt, size_t max_depth) {
    ASSERT_POINT_TYPE(T);
    return 0;
}

template class TreeNode<float, 2, 4>;               // 2D Quad-tree
template class TreeNode<float, 3, 8>;               // 3D Octree
template class TreeNode<double, 2, 4>;              // 2D Quad-tree
template class TreeNode<double, 3, 8>;              // 3D Octree

} // namespace name
