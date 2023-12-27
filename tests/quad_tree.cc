#include <iostream>
#include "utils/utils.h"
#include "trees/StaticMultiTree.h"

using namespace scds;

int main() {   
    StaticMultiTree<float, 2, 4> quad_tree(Point2f(0.5, 0.5), Point2f(0.5, 0.5), 32, 1);
    quad_tree.insert(Point2f(0.4359949 ,0.02592623));
    quad_tree.insert(Point2f(0.5496625 ,0.4353224));
    quad_tree.insert(Point2f(0.4203678 ,0.3303348));
    quad_tree.insert(Point2f(0.20464863,0.619271));
    quad_tree.insert(Point2f(0.29965466,0.2668273));
    quad_tree.insert(Point2f(0.6211338 ,0.5291421));
    quad_tree.insert(Point2f(0.13457994,0.5135781));

    std::vector<Point2f> nns, nns_bf;
    auto query = Point2f(0.23, 0.23);
    quad_tree.search_nn(query, nns, 8, 0.3);
    std::cout << "Tree size: " << quad_tree.size() << std::endl;
    std::cout << "Nearest neighbor number: "<< nns.size() << std::endl;
    for (const auto& pt: nns) {
        std::cout << pt << std::endl;
    }
    quad_tree.search_nn_bf(query, nns_bf, 8, 0.3);
    std::cout << "Nearest neighbor number (BF): "<< nns_bf.size() << std::endl;
    for (const auto& pt: nns_bf) {
        std::cout << pt << std::endl;
    }
    return 0;
}