# Minimal verions of cmake used
cmake_minimum_required (VERSION 3.0)

# The subproject
project (test_heat_fdtd)

# Include gtest library
find_package( GTest REQUIRED)
include_directories( ${GTEST_INCLUDE_DIR} )
enable_testing()
# TODO develope tests

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Test sources
set (SRC_TEST
    unit_test/test_fdtd_h1d.cpp
    )

# Executable test
add_executable( fdtd_test ${SRC_TEST} )
target_link_libraries( fdtd_test HeatFDTD)# ${GTEST_LIBRARIES} pthread )
set_target_properties( fdtd_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test/unit_test" )
