#include <iostream>
#include "utils/utils.h"
#include "trees/StaticMultiTree.h"

using namespace scds;

int main() {   
    StaticMultiTree<float, 2, 4> quad_tree(Point2f(0.5, 0.5), Point2f(0.5, 0.5), 16, 1);
    quad_tree.insert(Point2f(0.41, 0.32));
    quad_tree.insert(Point2f(0.2, 0.61));
    quad_tree.insert(Point2f(0.23, 0.31));
    quad_tree.insert(Point2f(0.43, 0.33));
    quad_tree.insert(Point2f(0.19, 0.13));
    quad_tree.insert(Point2f(0.67, 0.32));
    std::cout << quad_tree.size() << std::endl;
    return 0;
}