#pragma once
/**
 * Octree CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/
#include "Tree.h"

namespace scds {

/**
 * static octree for static scenes
*/

template <typename T>
class StaticOctree: public TreeBase<T> {
    
};

}       // end namespace scds