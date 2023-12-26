#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "trees/StaticMultiTree.h"

namespace py = pybind11;

PYBIND11_MODULE(smt, m) {

    m.doc() = "Static multi-tree python binding (Quad/Oct trees)\n";

    const char* SIZE_DOC = \
        "size():\n"
        "Obtaining the number of points currently stored in the tree"
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

    const char* SEARCH_NN_BF_DOC = \
        "search_nn_bf(pt: np.ndarray, k: int, radius: [float, double]):\n"
        "(Brute-force) Find at most k nearest neighbors with specified (see `search_nn` for more info)"
    ;

    py::class_<scds::StaticMultiTree<float, 2, 4>>(m, "QuadTreef")
        .def(py::init<const py::array_t<float>&, size_t, size_t>())
        .def(py::init<const py::array_t<float>&, size_t, float, size_t, size_t>())
        .def("size", &scds::StaticMultiTree<float, 2, 4>::size_py, SIZE_DOC)
        .def("insert", &scds::StaticMultiTree<float, 2, 4>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::StaticMultiTree<float, 2, 4>::search_nn_py, SEARCH_NN_DOC)
        .def("search_nn_bf", &scds::StaticMultiTree<float, 2, 4>::search_nn_bf_py, SEARCH_NN_BF_DOC);
    

    py::class_<scds::StaticMultiTree<float, 3, 8>>(m, "OctTreef")
        .def(py::init<const py::array_t<float>&, size_t, size_t>())
        .def(py::init<const py::array_t<float>&, size_t, float, size_t, size_t>())
        .def("size", &scds::StaticMultiTree<float, 3, 8>::size_py, SIZE_DOC)
        .def("insert", &scds::StaticMultiTree<float, 3, 8>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::StaticMultiTree<float, 3, 8>::search_nn_py, SEARCH_NN_DOC)
        .def("search_nn_bf", &scds::StaticMultiTree<float, 3, 8>::search_nn_bf_py, SEARCH_NN_BF_DOC);

    py::class_<scds::StaticMultiTree<double, 2, 4>>(m, "QuadTreed")
        .def(py::init<const py::array_t<double>&, size_t, size_t>())
        .def(py::init<const py::array_t<double>&, size_t, double, size_t, size_t>())
        .def("size", &scds::StaticMultiTree<double, 2, 4>::size_py, SIZE_DOC)
        .def("insert", &scds::StaticMultiTree<double, 2, 4>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::StaticMultiTree<double, 2, 4>::search_nn_py, SEARCH_NN_DOC)
        .def("search_nn_bf", &scds::StaticMultiTree<double, 2, 4>::search_nn_bf_py, SEARCH_NN_BF_DOC);

    py::class_<scds::StaticMultiTree<double, 3, 8>>(m, "OctTreed")
        .def(py::init<const py::array_t<double>&, size_t, size_t>())
        .def(py::init<const py::array_t<double>&, size_t, double, size_t, size_t>())
        .def("size", &scds::StaticMultiTree<double, 3, 8>::size_py, SIZE_DOC)
        .def("insert", &scds::StaticMultiTree<double, 3, 8>::insert_py, INSERT_DOC)
        .def("search_nn", &scds::StaticMultiTree<double, 3, 8>::search_nn_py, SEARCH_NN_DOC)
        .def("search_nn_bf", &scds::StaticMultiTree<double, 3, 8>::search_nn_bf_py, SEARCH_NN_BF_DOC);
}

