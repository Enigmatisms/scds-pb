cmake_minimum_required(VERSION 3.10)
project(SCDS)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
IF(CMAKE_BUILD_TYPE EQUAL Release)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3 -march=native")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")
ELSE()
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check C++14 support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++17 COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
   add_definitions(-DCOMPILEDWITHC17)
   message(STATUS "Using flag -std=c++17.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Eigen3 REQUIRED)
if(NOT Eigen3_FOUND)
   message(FATAL_ERROR "Eigen3 not found.")
endif()

include_directories(
    .
    ./include
    ${EIGEN3_INCLUDE_DIR}
    /usr/local/include/
)

link_directories(/usr/local/lib)

add_subdirectory(app)
add_subdirectory(tests)