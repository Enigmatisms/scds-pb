#include <array>
#include <vector>
#include "Point.h"
namespace scds {

template<typename T, size_t Ndim, size_t Nchild>
class TreeNode {
public:

private:
    std::array<TreeNode*, Nchild> childs;
    Point<T, Ndim> center;
    Point<T, Ndim> size;
    bool is_leaf;
    int  num_point;
};

template<typename T>
using QdNode = TreeNode<T, 2, 4>;
template<typename T>
using OcNode = TreeNode<T, 3, 8>;

} // end namespace scds