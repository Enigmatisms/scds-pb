#include <iostream>
#include "utils/TreeNode.h"
#include "utils/utils.h"

using namespace scds;

int main()
{   
    Point3f tl{0, 0, 0};
    Point3f br{1, 2, 3};
    std::vector<Point3f> memory_arena;
    std::unordered_set<size_t> subset_idxs;
    OcNode<float> node{tl, br, std::move(subset_idxs), nullptr, std::make_shared<std::vector<Point3f>>(memory_arena)};
    return 0;
}