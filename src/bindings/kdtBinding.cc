#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "trees/KDTree.h"

namespace py = pybind11;

PYBIND11_MODULE(kdt, kdtree_module) {

    kdtree_module.doc() = "k-d tree python binding\n";

    const char* SIZE_DOC = \
        "size():\n"
        "Obtaining the number of points currently stored in the tree"
    ;

    const char* DEPTH_DOC = \
        "depth():\n"
        "Obtaining the current depth of the tree"
    ;

    const char* INSERT_DOC = \
        "insert(pt: np.ndarray):\n"
        "Insert point(s) to the tree\n"
        "pt: point(s) to insert. Should be of shape (Ndim) or (N_points, Ndim)"
    ;

    const char* SEARCH_NN_DOC = \
        "search_nn(pt: np.ndarray, k: int, radius: [float, double]):\n"
        "Find at most k nearest neighbors with specified\n"
        "pt: point around which the search is performed. Should be of shape (Ndim) or (1, Ndim)\n"
        "k: number of nearest neighbors to keep (maximum, can be fewer than this if there aren't enough neighbors)\n"
        "radius: spatial tree search radius\n\n"
        "return: (N_found, Ndim) np.ndarray with N_found neighbors, note that\n"
        "the points returned are ordered by distance to `pt` (descending)"
    ;

    const char* BUILD_TREE_DOC = \
        "build_tree(pt: np.ndarray):\n"
        "Build the k-d tree all at once with a large batch of points\n"
        "pt: point around which the search is performed. Should be of shape (Ndim) or (1, Ndim)\n"
    ;

    const char* TREE_STRUCTURE_DOC = \
        "search_nn(pt: np.ndarray, k: int, radius: [float, double]):\n"
        "Find at most k nearest neighbors with specified\n"
        "pt: point around which the search is performed. Should be of shape (Ndim) or (1, Ndim)\n"
        "k: number of nearest neighbors to keep (maximum, can be fewer than this if there aren't enough neighbors)\n"
        "radius: spatial tree search radius\n\n"
        "return: (N_found, Ndim) np.ndarray with N_found neighbors, note that\n"
        "the points returned are ordered by distance to `pt` (descending)"
    ;

    py::class_<scds::KDTree<float, 2>>(kdtree_module, "KDTree2f")
        .def(py::init<const py::array_t<float>&, int, int, int, float>(), py::arg("bbox_info"), py::arg("max_depth") = 0, py::arg("node_max_point_num") = 1, py::arg("k") = 1, py::arg("radius") = 0)
        .def(py::init<const py::array_t<float>&, size_t, float, int, int, int, float>())
        .def("size", &scds::KDTree<float, 2>::size_py, SIZE_DOC)
        .def("depth", &scds::KDTree<float, 2>::depth_py, DEPTH_DOC)
        .def("insert", &scds::KDTree<float, 2>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::KDTree<float, 2>::search_nn_py, SEARCH_NN_DOC)
        .def("build_tree", &scds::KDTree<float, 2>::build_tree_py, BUILD_TREE_DOC)
        .def("tree_structure", &scds::KDTree<float, 2>::get_tree_structure, TREE_STRUCTURE_DOC);
    

    py::class_<scds::KDTree<float, 3>>(kdtree_module, "KDTree3f")
        .def(py::init<const py::array_t<float>&, int, int, int, float>())
        .def(py::init<const py::array_t<float>&, size_t, float, int, int, int, float>())
        .def("size", &scds::KDTree<float, 3>::size_py, SIZE_DOC)
        .def("depth", &scds::KDTree<float, 3>::depth_py, DEPTH_DOC)
        .def("insert", &scds::KDTree<float, 3>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::KDTree<float, 3>::search_nn_py, SEARCH_NN_DOC)
        .def("build_tree", &scds::KDTree<float, 3>::build_tree_py, BUILD_TREE_DOC)
        .def("tree_structure", &scds::KDTree<float, 3>::get_tree_structure, TREE_STRUCTURE_DOC);

    py::class_<scds::KDTree<double, 2>>(kdtree_module, "KDTree2d")
        .def(py::init<const py::array_t<double>&, int, int, int, double>())
        .def(py::init<const py::array_t<double>&, size_t, double, int, int, int, float>())
        .def("size", &scds::KDTree<double, 2>::size_py, SIZE_DOC)
        .def("depth", &scds::KDTree<double, 2>::depth_py, DEPTH_DOC)
        .def("insert", &scds::KDTree<double, 2>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::KDTree<double, 2>::search_nn_py, SEARCH_NN_DOC)
        .def("build_tree", &scds::KDTree<double, 2>::build_tree_py, BUILD_TREE_DOC)
        .def("tree_structure", &scds::KDTree<double, 2>::get_tree_structure, TREE_STRUCTURE_DOC);

    py::class_<scds::KDTree<double, 3>>(kdtree_module, "KDTree3d")
        .def(py::init<const py::array_t<double>&, int, int, int, double>())
        .def(py::init<const py::array_t<double>&, size_t, double, int, int, int, float>())
        .def("size", &scds::KDTree<double, 3>::size_py, SIZE_DOC)
        .def("depth", &scds::KDTree<double, 3>::depth_py, DEPTH_DOC)
        .def("insert", &scds::KDTree<double, 3>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::KDTree<double, 3>::search_nn_py, SEARCH_NN_DOC)
        .def("build_tree", &scds::KDTree<double, 3>::build_tree_py, BUILD_TREE_DOC)
        .def("tree_structure", &scds::KDTree<double, 3>::get_tree_structure, TREE_STRUCTURE_DOC);
}

