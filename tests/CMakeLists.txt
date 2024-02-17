# staticMultiTree (SMT) python binding


# ============== pybind11 python module declaration ===================
set(SMT_TEST_SOURCES 
    ../src/trees/StaticMultiTree.cc
    ../src/bindings/smtBinding.cc
    ../src/utils/TreeNode.cc
    ../src/utils/stats.cc
)

pybind11_add_module(smt ${SMT_TEST_SOURCES})
target_include_directories(smt PRIVATE 
    ../ext/fmt/include/
    ${Python_INCLUDE_DIRS}
)
target_link_libraries(smt PRIVATE fmt::fmt pybind11::module ${PYTHON_LIBRARIES} Python::Python)
set_target_properties(smt PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")

set(KDT_TEST_SOURCES 
    ../src/trees/KDTree.cc
    ../src/bindings/kdtBinding.cc
    ../src/utils/BinTreeNode.cc
    ../src/utils/stats.cc
)

pybind11_add_module(kdt ${KDT_TEST_SOURCES})
target_include_directories(kdt PRIVATE 
    ../ext/fmt/include/
    ${Python_INCLUDE_DIRS}
)
target_link_libraries(kdt PRIVATE fmt::fmt pybind11::module ${PYTHON_LIBRARIES} Python::Python)
set_target_properties(kdt PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")
# ============== end of module declaration ===================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

aux_source_directory(../src/utils UTILS_SRC)
aux_source_directory(../src/trees TREES_SRC)
add_executable(point_test point_test.cc)

# quad tree unit test
add_executable(tree_node tree_node.cc
    ${UTILS_SRC}
)

# quad tree unit test
add_executable(quad_tree quad_tree.cc
    ${UTILS_SRC} ${TREES_SRC}
)

add_executable(kdtree kdtree.cc
    ${UTILS_SRC} ${TREES_SRC}
)

target_link_libraries(quad_tree PRIVATE fmt::fmt pybind11::module ${PYTHON_LIBRARIES} Python::Python)
target_link_libraries(kdtree PRIVATE fmt::fmt pybind11::module ${PYTHON_LIBRARIES} Python::Python)

