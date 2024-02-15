#pragma once
/**
 * spatial tree (regular parition) CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/

#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "utils/BinTreeNode.h"
#include "utils/stats.h"

namespace scds {

#define KDT_MEMORY_PROFILE

#ifdef KDT_MEMORY_PROFILE

STAT_MEMORY_COUNTER("KDT/Total points", pointBytes);

#endif //KDT_MEMORY_PROFILE

// comparing functor
template<typename T>
struct DistanceComp {
    bool operator()(const std::pair<int, T>& pr1, const std::pair<int, T>& pr2) const {
        return pr1.second < pr2.second;
    }
};

/**
 * @brief KD-tree efficient nearest neighbor earch structure
*/
template <typename T, size_t Ndim>
class KDTree {

public:
    using Node     = BinTreeNode<T, Ndim>;
    using Pointx   = Point<T, Ndim>;
    using PointVec = std::vector<Point<T, Ndim>>;
    using QueueType = std::priority_queue<std::pair<int, T>, std::vector<std::pair<int, T>>, DistanceComp<T>>;
    // pybind initializer (1)
    KDTree(const pybind11::array_t<T>& bbox_info, size_t max_depth = 0, size_t node_max_point_num = 0, int k = 1, T radius = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM)),
        k(k), radius(radius), max_heap(DistanceComp<T>{})
    {
        all_pts->reserve(64);
        root = std::make_shared<Node>(
            Pointx::from_pointer(bbox_info.data()), 
            Pointx::from_pointer(bbox_info.data() + Ndim),
            std::weak_ptr<Node>()
        );
    }
    
    template <typename Ptype1, typename Ptype2>
    KDTree(Ptype1&& center, Ptype2&& half_size, size_t max_depth = 0, size_t node_max_point_num = 0, int k = 1, T radius = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM)),
        k(k), radius(radius), max_heap(DistanceComp<T>{})
    {
        all_pts->reserve(64);
        root = std::make_shared<Node>(
            std::forward<Ptype1>(center), 
            std::forward<Ptype2>(half_size),
            std::weak_ptr<Node>()
        );
    }

    // pybind initializer (2)
    KDTree(
        const pybind11::array_t<T>& points, size_t num_points, T border = 0, 
        size_t max_depth = 0, size_t node_max_point_num = 0, int k = 1, T radius = 0
    ):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM)),
        k(k), radius(radius), max_heap(DistanceComp<T>{})
    {
        all_pts->reserve(num_points);
        const T* ptr = points.data();
        Pointx min_range = Pointx::from_pointer(ptr);
        Pointx max_range = min_range;
        for (size_t point_cnt = 0; point_cnt < num_points; ptr += Ndim, point_cnt++) {
            auto pt = Pointx::from_pointer(ptr);
            max_range.max_inplace(pt);
            min_range.min_inplace(pt);
        }
        max_range = (max_range + min_range) / 2;    // center
        min_range = max_range - min_range + static_cast<float>(border);          // half_size

        root = std::make_shared<Node>(
            std::move(max_range), 
            std::move(min_range),
            std::weak_ptr<Node>()
        );

        ptr = points.data();
        for (size_t point_cnt = 0; point_cnt < num_points; ptr += Ndim, point_cnt++) {
            auto pt = Pointx::from_pointer(ptr);
            insert(pt);
        }
    }

    template <typename PointVecType>
    KDTree(PointVecType&& points, T border = 0, size_t max_depth = 0, size_t node_max_point_num = 0, int k = 1, T radius = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM)),
        k(k), radius(radius), max_heap(DistanceComp<T>{})
    {
        all_pts->reserve(points.size());
        Pointx min_range = points.front();
        Pointx max_range = points.front();
        for (const auto& pt: points) {
            max_range.max_inplace(pt);
            min_range.min_inplace(pt);
        }
        max_range = (max_range + min_range) / 2;    // center
        min_range = max_range - min_range + static_cast<float>(border);          // half_size

        root = std::make_shared<Node>(
            std::move(max_range), 
            std::move(min_range),
            std::weak_ptr<Node>()
        );

        for (const auto& pt: points)
            insert(pt);
    }
public:
    // insert point in the tree
    void insert(const Pointx& pt);

    /**
     * @brief Nearest neighbor search in the K-D-tree, note that k-d tree does not rely on search radius strictly
     * @param pt:     point to search around
     * @param nn:    resulting points are stored here
     * 
     * Actually, quad/oct tree is not suitable for KNN (specifying K), since tree will not be easy to prune
     * so here, radius will be a compulsory parameter
    */
    void search_nn(const Pointx& pt, PointVec& nn, int k = 1, T radius = 0) const;

    void recursive_solve(const Pointx& pt, QueueType& queue, std::shared_ptr<Node> cur_node) const;

    int size() const {
        return root ? root->num_points : 0;
    }

    int depth() const {
        return this->tree_depth;
    }
public:     // python binding
    void insert_py(const pybind11::array_t<T>& pt);

    pybind11::array_t<T> search_nn_py(const pybind11::array_t<T>& pt, int k = 1, T radius = 0) const;

    // get the structure of the tree
    pybind11::tuple get_tree_structure() const;

    int size_py() const {return static_cast<int>(root ? root->num_points : 0);}
    int depth_py() const {return static_cast<int>(this->tree_depth);}
protected:
    static constexpr size_t valid_num_check(size_t value, size_t max_num) {
        value = (value == 0) ? max_num : value;
        return std::min(max_num, value);
    }
protected:
    std::shared_ptr<PointVec> all_pts;
    std::shared_ptr<Node> root;
private:
    size_t tree_depth{0};
    // maximum depth of the tree (1st priority)
    const size_t max_depth;
    // maximum number of point in a node (2nd priority)
    const size_t node_max_point_num;
    // number of nearest points to extract
    int k;
    // radius from which the nearest points are extract
    T radius;

    static constexpr size_t MAX_DEPTH    = 32;                         // physical barrier
    static constexpr size_t MAX_NODE_NUM = 64;                         // physical barrier
};

/**
 * @todo: ...

*/
template <typename T>
using KDTree2 = KDTree<T, 2>;

/**
 * @todo: ...
*/
template <typename T>
using KDTree3 = KDTree<T, 3>;

}       // end namespace scds