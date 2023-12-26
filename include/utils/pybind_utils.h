#pragma once
#include "utils.h"
#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace scds {

#ifndef NO_PYBIND_SCOPE
#ifdef SCDS_DEBUG
template <typename T>
bool pyArrayShapeCheck(const pybind11::array_t<T>& pt, int Ndim) {
    int ndim = pt.ndim();
    if (ndim == 1) return true;
    else {
        if (pt.shape()[1] == 1) return true;
        return false;
    }
}

#define SMT_ARRAY_DTYPE_CHECK(...)

#else   // SCDS_DEBUG
template <typename T>
bool pyArrayShapeCheck(const pybind11::array_t<T>& pt, int Ndim) {
    int ndim = pt.ndim();
    if (ndim > 2)
        SCDS_RUNTIME_ERROR("Ndim of the input array is greater than 2, currently: {}", ndim);
    else if (ndim == 1) {
        if (pt.shape()[0] != Ndim)
            SCDS_RUNTIME_ERROR("Point dimension miss match: Ndim = {0}, shape: {1}", Ndim, pt.shape()[0]);
        return true;
    } else {
        if (pt.shape()[1] == 1) {
            if (pt.shape()[0] != Ndim)
                SCDS_RUNTIME_ERROR("Point dimension miss match: Ndim = {0}, shape: {1}", Ndim, pt.shape()[0]);
            return true;
        } else {
            if (pt.shape()[1] != Ndim)
                SCDS_RUNTIME_ERROR("Point dimension miss match: Ndim = {0}, shape: {1}", Ndim, pt.shape()[1]);
            return false;
        }
    }
}

// 
#define SMT_ARRAY_DTYPE_CHECK(pt, Type) \
if (pt.dtype().is(pybind11::dtype::of<Type>())) { \
    if constexpr (std::is_same_v<Type, float>) \
        SCDS_RUNTIME_ERROR("Input array must have dtype float32."); \
    else \
        SCDS_RUNTIME_ERROR("Input array must have dtype float64."); \
}
#endif // SCDS_DEBUG

template <typename T>
pybind11::array_t<T> create_array_2d(int row, int col) {
    return pybind11::array_t<T>({row, col}, {col * sizeof(T), sizeof(T)});
}

template <typename T>
pybind11::array_t<T> create_array_1d(int col) {
    return pybind11::array_t<T>({col}, {sizeof(T)});
}

#endif // NO_PYBIND
    
} // namespace scds
