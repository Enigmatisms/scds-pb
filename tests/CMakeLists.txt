aux_source_directory(
    ../src SRC
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

add_executable(point_test point_test.cc
    ${SRC}
)
