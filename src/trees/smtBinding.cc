#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include "trees/StaticMultiTree.h"

namespace py = pybind11;

PYBIND11_MODULE(smt, m) {

    m.doc() = "Static multi-tree python binding (Quad/Oct trees)\n";

    py::class_<scds::StaticMultiTree<float, 2, 4>>(m, "QuadTreef")
        .def(py::init<const py::array_t<float>&, const py::array_t<float>&, size_t, size_t>())
        .def(py::init<const py::array_t<float>&, float, size_t, size_t>())
        .def("insert", &scds::StaticMultiTree<float, 2, 4>::insert)
        .def("search_nn", &scds::StaticMultiTree<float, 2, 4>::search_nn)
        .def("search_nn_bf", &scds::StaticMultiTree<float, 2, 4>::search_nn_bf);
    

    py::class_<scds::StaticMultiTree<float, 3, 8>>(m, "OctTreef")
        .def(py::init<const py::array_t<float>&, const py::array_t<float>&, size_t, size_t>())
        .def(py::init<const py::array_t<float>&, float, size_t, size_t>())
        .def("insert", &scds::StaticMultiTree<float, 3, 8>::insert)
        .def("search_nn", &scds::StaticMultiTree<float, 3, 8>::search_nn)
        .def("search_nn_bf", &scds::StaticMultiTree<float, 3, 8>::search_nn_bf);

    py::class_<scds::StaticMultiTree<double, 2, 4>>(m, "QuadTreef")
        .def(py::init<const py::array_t<double>&, const py::array_t<double>&, size_t, size_t>())
        .def(py::init<const py::array_t<double>&, double, size_t, size_t>())
        .def("insert", &scds::StaticMultiTree<double, 2, 4>::insert)
        .def("search_nn", &scds::StaticMultiTree<double, 2, 4>::search_nn)
        .def("search_nn_bf", &scds::StaticMultiTree<double, 2, 4>::search_nn_bf);

    py::class_<scds::StaticMultiTree<double, 3, 8>>(m, "OctTreef")
        .def(py::init<const py::array_t<double>&, const py::array_t<double>&, size_t, size_t>())
        .def(py::init<const py::array_t<double>&, double, size_t, size_t>())
        .def("insert", &scds::StaticMultiTree<double, 3, 8>::insert)
        .def("search_nn", &scds::StaticMultiTree<double, 3, 8>::search_nn)
        .def("search_nn_bf", &scds::StaticMultiTree<double, 3, 8>::search_nn_bf);
}

