#include "utils/BinTreeNode.h"

namespace scds {

template<typename Ty, size_t Ndim>
SplitAxis BinTreeNode<Ty, Ndim>::max_extent_axis(const std::vector<Pointx>& pts) const {
    Pointx maxi = pts[sub_idxs->front()], mini = maxi;
    for (size_t i = 1; i < sub_idxs->size(); i++) {
        const Pointx& pt = pts[sub_idxs->at(i)];
        maxi = maxi.max(pt);
        mini = mini.min(pt);
    }
    Pointx diff = maxi - mini;
    int max_axis  = 0;
    Ty max_extent = diff[0]; 
    for (int i = 1; i < 3; i++) {
        if (diff[i] > max_extent) {
            max_extent = diff[i];
            max_axis   = i;
        } 
    }
    return SplitAxis(max_axis);
}

template<typename Ty, size_t Ndim>
void BinTreeNode<Ty, Ndim>::split_leaf_node(const std::vector<Pointx>& pts) {
    // get split axis
    SplitAxis axis   = max_extent_axis(pts);
    this->split_axis = axis; 
    // get split pos using nth_element
    auto center_it = pts.begin() + (pts.size() >> 1);
    auto comp_op = [axis](const Pointx& p1, const Pointx& p2) {return p1[axis] < p2[axis];};
    // this is not correct (huge problem)
    std::vector<Pointx>::iterator elem_1 = std::nth_element(pts.begin(), center_it, pts.end(), comp_op),
         elem_2 = std::nth_element(pts.begin(), center_it + 1, pts.end(), comp_op);
    this->split_pos = 0.5 * ((*elem_1)[axis] + (*elem_2)[axis]);
    std::vector<int> lcontainer, rcontainer;
    for (int idx: *sub_idxs) {
        if (pts[idx][axis] < this->split_pos)
            lcontainer.push_back(idx);
        else
            rcontainer.push_back(idx);
    }
    Ty lsize = (half_size[axis] - center[axis] + split_pos) / 2, rsize = half_size[axis] - lsize,
       lpos  = split_pos - lsize, rpos = split_pos + rsize;
    Pointx l_ctr = this->center, r_ctr = this->center,
           l_size = this->half_size, r_size = this->half_size;
    l_ctr[axis]  = lpos;
    l_size[axis] = lsize;
    r_ctr[axis]  = rpos;
    r_size[axis] = rsize;
    add_child(std::move(l_ctr), std::move(l_size), std::move(lcontainer), true);
    add_child(std::move(r_ctr), std::move(r_size), std::move(rcontainer), false);
}

template<typename Ty, size_t Ndim>
size_t BinTreeNode<Ty, Ndim>::get_size() const {
    size_t size = sizeof(*this);  // Size of the current instance
    if (is_leaf())
        size += sub_idxs->size() * sizeof(int);
    size += sizeof(_lchild);
    size += sizeof(_rchild);
    return size;
}

template class BinTreeNode<float, 2>;               // 2D KD-tree (vector container)
template class BinTreeNode<float, 3>;               // 3D KD-tree (vector container)
template class BinTreeNode<double, 2>;              // 2D KD-tree (vector container)
template class BinTreeNode<double, 3>;              // 3D KD-tree (vector container)

} // namespace name
