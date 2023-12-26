#include <iostream>
#include "utils/utils.h"
#include "utils/TreeNode.h"

#define NO_PYBIND_SCOPE

using namespace scds;

int main() {   
    std::shared_ptr<TreeNode<float, 2, 4>> root_node;
    root_node = std::make_shared<TreeNode<float, 2, 4>>(Point2f(0.5, 0.5), Point2f(0.5, 0.5), std::weak_ptr<TreeNode<float, 2, 4>>(), std::make_shared<std::vector<Point2f>>());
    auto child = root_node->try_get_child(0);
    std::cout << child->num_points() << std::endl;
    return 0;
}