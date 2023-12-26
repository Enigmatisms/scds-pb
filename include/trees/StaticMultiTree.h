#pragma once
/**
 * spatial tree (regular parition) CPU / GPU implementation
 * @author: Qianyue He
 * @date:   2023-12-14
*/

#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "utils/TreeNode.h"

namespace scds {

/**
 * @brief Static spatial tree (quad / oct) for static scenes.
*/
template <typename T, size_t Ndim, size_t Nchild>
class StaticMultiTree {
public:
    using Node     = TreeNode<T, Ndim, Nchild>;
    using Pointx   = Point<T, Ndim>;
    using PointVec = std::vector<Point<T, Ndim>>;

    // pybind initializer (1)
    StaticMultiTree(const pybind11::array_t<T>& bbox_info, size_t max_depth = 0, size_t node_max_point_num = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM))
    {
        all_pts->reserve(64);
        root = std::make_shared<Node>(
            Pointx::from_pointer(bbox_info.data()),
            Pointx::from_pointer(bbox_info.data() + Ndim),
            std::weak_ptr<Node>(), 
            all_pts
        );
    }
    
    template <typename Ptype1, typename Ptype2>
    StaticMultiTree(Ptype1&& center, Ptype2&& half_size, size_t max_depth = 0, size_t node_max_point_num = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM))
    {
        all_pts->reserve(64);
        root = std::make_shared<Node>(
            std::forward<Ptype1>(center), 
            std::forward<Ptype2>(half_size),
            std::weak_ptr<Node>(), 
            all_pts
        );
    }

    // pybind initializer (2)
    StaticMultiTree(const pybind11::array_t<T>& points, size_t num_points, T border = 0, size_t max_depth = 0, size_t node_max_point_num = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM))
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
            std::weak_ptr<Node>(),
            all_pts
        );

        ptr = points.data();
        for (size_t point_cnt = 0; point_cnt < num_points; ptr += Ndim, point_cnt++) {
            auto pt = Pointx::from_pointer(ptr);
            insert(pt);
        }
    }

    template <typename PointVecType>
    StaticMultiTree(PointVecType&& points, T border = 0, size_t max_depth = 0, size_t node_max_point_num = 0):
        all_pts(std::make_shared<PointVec>()),
        max_depth(valid_num_check(max_depth, MAX_DEPTH)), 
        node_max_point_num(valid_num_check(node_max_point_num, MAX_NODE_NUM))
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
            std::weak_ptr<Node>(),
            all_pts
        );

        for (const auto& pt: points)
            insert(pt);
    }
public:
    // insert point in the tree
    void insert(const Pointx& pt);

    /**
     * @brief Nearest neighbor search in the tree
     * @param pt:     point to search around
     * @param knn:    resulting points are stored here
     * @param k:      number of nearest points
     * @param radius: we can specify the radius to keep the point in this radius. 
     * When radius search is used, be sure to set k = 0
     * 
     * Actually, quad/oct tree is not suitable for KNN (specifying K), since tree will not be easy to prune
     * so here, radius will be a compulsory parameter
    */
    void search_nn(const Pointx& pt, PointVec& nn, size_t k = 1, T radius = static_cast<T>(0)) const;

    // same utility but implemented via brute force searching
    void search_nn_bf(const Pointx& pt, PointVec& nn, size_t k = 1, T radius = static_cast<T>(0)) const;

    size_t size() const {
        return root ? root->num_points() : 0;
    }
public:     // python binding
    void insert_py(const pybind11::array_t<T>& pt);

    pybind11::array_t<T> search_nn_py(const pybind11::array_t<T>& pt, int k = 1, T radius = static_cast<T>(0)) const;

    pybind11::array_t<T> search_nn_bf_py(const pybind11::array_t<T>& pt, int k = 1, T radius = static_cast<T>(0)) const;

    int size_py() const {return static_cast<int>(root ? root->num_points() : 0);}
protected:
    /**
     * @brief decide which child the node is in (via bit operation), for example, in 2D:
     * @note    y
     * @note 01 | 11
     * @note ---.--- x
     * @note 00 | 10
    */
    template <typename PointType>
    static size_t which_child(std::shared_ptr<Node> node, PointType&& pt) {
        STATIC_ASSERT_POINT_TYPE(PointType);
        auto judge = ((pt - node->center) > 0).to_u64();
        return (judge[0] << 1) + judge[1];
    }

    static constexpr size_t valid_num_check(size_t value, size_t max_num) {
        value = (value == 0) ? max_num : value;
        return std::min(max_num, value);
    }
protected:
    std::shared_ptr<PointVec> all_pts;
    std::shared_ptr<Node> root;
private:
    // maximum depth of the tree (1st priority)
    const size_t max_depth;
    // maximum number of point in a node (2nd priority)
    const size_t node_max_point_num;

    static constexpr size_t MAX_DEPTH    = 32;                         // physical barrier
    static constexpr size_t MAX_NODE_NUM = 64;                         // physical barrier
};

/**
 * QuadTree (2D). Initialize the tree by passing:
 * center: center of the bounding box
 * half_size: half size of the bounding box
 * max_depth: 0 by default, which means no maximum depth limit
 * node_max_point_num: maximum number of point in a leaf node
*/
template <typename T>
using QuadTree = StaticMultiTree<T, 2, 4>;

/**
 * OctTree (3D). Initialize the tree by passing:
 * center: center of the bounding box
 * half_size: half size of the bounding box
 * max_depth: 0 by default, which means no maximum depth limit
 * node_max_point_num: maximum number of point in a leaf node
*/
template <typename T>
using OctTree = StaticMultiTree<T, 3, 8>;

}       // end namespace scds