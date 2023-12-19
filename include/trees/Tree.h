#pragma once
/**
 * Tree base class implementation
 * 
*/

#include "utils/Point.h"

namespace scds {

template<typename T>
class TreeBase {
public:
    TreeBase();
public:
    template <typename PointType>
    virtual void insert(PointType&& pt) = 0;

    virtual void tree_builder() = 0;
};

}   // end of namespace scds