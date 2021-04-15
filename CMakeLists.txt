cmake_minimum_required(VERSION 3.5)
project(fitpackpp VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
option(BUILD_EXAMPLES "Build examples" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(compile_options)

add_subdirectory(fitpack)
add_subdirectory(fitpackpp)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
