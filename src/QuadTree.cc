#include "trees/Quadtree.h"

namespace scds {

template<typename Ty>
template<typename PointType>
void StaticQuadTree<Ty>::insert(PointType&& pt, size_t max_depth) {
    ASSERT_POINT_TYPE(PointType);
    // 'recursive-like' tree-build
    return 0;
}


}   // end namespace scds