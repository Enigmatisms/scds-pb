# staticMultiTree (SMT) python binding

aux_source_directory(../src/utils UTILS_SRC)
aux_source_directory(../src/trees TREES_SRC)

# quad tree unit test
add_executable(benchmark benchmark.cc
    ${UTILS_SRC} ${TREES_SRC}
)

target_link_libraries(benchmark PRIVATE fmt::fmt pybind11::module ${PYTHON_LIBRARIES} Python::Python)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})