#include <iostream>
#include <random>
#include "utils/stats.h"
#include "utils/utils.h"
#include "trees/StaticMultiTree.h"

using namespace scds;

int main(int argc, char** argv) {   
    InitProfiler();
    StaticMultiTree<float, 2, 4> quad_tree(Point2f(0.5, 0.5), Point2f(0.5, 0.5), 32, 1);

    size_t num_points = 2000000;
    if (argc > 1) {
        num_points = static_cast<size_t>(atoi(argv[1]));
        num_points = std::max(8ul, num_points);
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0, 1.0);
    for (size_t i = 0; i < num_points; i++) {
        float rand_x = dis(gen);
        float rand_y = dis(gen);
        quad_tree.insert(Point2f(rand_x , rand_y));
    }
    // ============== smaller scale test ===============
    // quad_tree.insert(Point2f(0.4359949 ,0.02592623));
    // quad_tree.insert(Point2f(0.5496625 ,0.4353224));
    // quad_tree.insert(Point2f(0.4203678 ,0.3303348));
    // quad_tree.insert(Point2f(0.20464863,0.619271));
    // quad_tree.insert(Point2f(0.29965466,0.2668273));
    // quad_tree.insert(Point2f(0.6211338 ,0.5291421));
    // quad_tree.insert(Point2f(0.13457994,0.5135781));

    std::vector<Point2f> nns, nns_bf;
    auto query = Point2f(0.23, 0.23);
    quad_tree.search_nn(query, nns, 8, 0.05);
    std::cout << "Tree size: " << quad_tree.size() << std::endl;
    std::cout << "Nearest neighbor number: "<< nns.size() << std::endl;
    
    quad_tree.search_nn_bf(query, nns_bf, 8, 0.05);
    std::cout << "Nearest neighbor number (BF): "<< nns_bf.size() << std::endl;
    std::cout << "Tree depth before destruct: " << quad_tree.depth() << std::endl;

    ReportThreadStats();    
    PrintStats(stdout);
    ReportProfilerResults(stdout);

    ClearStats();
    ClearProfiler();
    CleanupProfiler();
    return 0;
}