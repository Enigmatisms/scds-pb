set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

aux_source_directory(../src/utils UTILS_SRC)
aux_source_directory(../src/trees TREES_SRC)
add_executable(point_test point_test.cc)

add_executable(tree_node tree_node.cc
    ${UTILS_SRC}
)

add_executable(quad_tree quad_tree.cc
    ${UTILS_SRC} ${TREES_SRC}
)
