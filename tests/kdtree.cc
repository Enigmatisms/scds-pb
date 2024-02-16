#include <iostream>
#include <random>
#include "utils/stats.h"
#include "utils/utils.h"
#include "trees/KDTree.h"

using namespace scds;

int main(int argc, char** argv) {   
    InitProfiler();
    KDTree2<float> kdtree(Point2f(0.5, 0.5), Point2f(0.5, 0.5), 32, 4);

    // size_t num_points = 2000000;
    // if (argc > 1) {
    //     num_points = static_cast<size_t>(atoi(argv[1]));
    //     num_points = std::max(8ul, num_points);
    // }
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<float> dis(0.0, 1.0);
    // for (size_t i = 0; i < num_points; i++) {
    //     float rand_x = dis(gen);
    //     float rand_y = dis(gen);
    //     kdtree.insert(Point2f(rand_x , rand_y));
    // }
    
    // std::vector<Point2f> nns, nns_bf;
    // auto query = Point2f(0.23, 0.23);
    // kdtree.search_nn(query, nns, 8, 0.05);
    // std::cout << "Tree size: " << kdtree.size() << std::endl;
    // std::cout << "Nearest neighbor number: "<< nns.size() << std::endl;
    
    // ReportThreadStats();    
    // PrintStats(stdout);
    // ReportProfilerResults(stdout);

    ClearStats();
    ClearProfiler();
    CleanupProfiler();
    return 0;
}