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
    for (size_t i = 1; i < Ndim; i++) {
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
    auto comp_op = [axis, &pts](int index_1, int index_2) {return pts[index_1][axis] < pts[index_2][axis];};

    int half_num = static_cast<int>(sub_idxs->size() - 1) >> 1;
    std::nth_element(sub_idxs->begin(), sub_idxs->begin() + half_num, sub_idxs->end(), comp_op);
    int index_1 = *(sub_idxs->begin() + half_num);
    std::nth_element(sub_idxs->begin(), sub_idxs->begin() + half_num + 1, sub_idxs->end(), comp_op);
    int index_2 = *(sub_idxs->begin() + half_num + 1);

    this->split_pos = (pts[index_1][axis] + pts[index_2][axis]) * static_cast<Ty>(0.5);

    std::vector<int> lcontainer, rcontainer;
    for (int idx: *sub_idxs) {
        if (pts[idx][axis] < this->split_pos)
            lcontainer.push_back(idx);
        else
            rcontainer.push_back(idx);
    }
    sub_idxs.reset(nullptr);
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
