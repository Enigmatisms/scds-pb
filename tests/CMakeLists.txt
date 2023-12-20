aux_source_directory(../src/utils UTILS_SRC)
aux_source_directory(../src/trees TREES_SRC)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

add_executable(point_test point_test.cc
    ${UTILS_SRC}
)

add_executable(tree_node tree_node.cc
    ${UTILS_SRC} ${TREES_SRC}
)
