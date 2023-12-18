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
    virtual void insert();
    virtual void tree_builder();
};

}   // end of namespace scds