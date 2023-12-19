aux_source_directory(
    ../src SRC
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(dummy dummy.cc
    ${SRC}
)
