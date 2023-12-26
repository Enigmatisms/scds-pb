#include <iostream>
#include "utils/utils.h"
#include "trees/StaticMultiTree.h"

using namespace scds;

int main() {   
    StaticMultiTree<float, 2, 4> quad_tree(Point2f(0.5, 0.5), Point2f(0.5, 0.5), 16, 1);
    quad_tree.insert(Point2f(0.2, 0.2));
    quad_tree.insert(Point2f(0.6, 0.6));
    quad_tree.insert(Point2f(0.2, 0.6));
    quad_tree.insert(Point2f(0.6, 0.2));
    quad_tree.insert(Point2f(0.4, 0.4));
    std::cout << quad_tree.size() << std::endl;
    return 0;
}