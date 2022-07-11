# Minimal version of cmake used
cmake_minimum_required (VERSION 3.0)

# The project
project (HeatFDTD)

# Set directories
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/src)

# Include directories
include_directories(src)
include_directories(${CMAKE_BINARY_DIR}/../Vector_Matrix/src)
include_directories(${CMAKE_BINARY_DIR}/../MaterialManagement/src)
include_directories(${CMAKE_BINARY_DIR}/../MaterialManagement/ThirdParty)

# Include library path
link_directories(${CMAKE_BINARY_DIR}/../MaterialManagement/bin)

# Compilation flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O3")

set(SRC
    src/heat_fdtd_1d.cpp
    )

set(HEADER
    src/heat_fdtd_1d.hpp
    )

add_library(HeatFDTD SHARED ${SRC} ${HEADER})
target_link_libraries( HeatFDTD MaterialManagement pthread )
#set_target_properties( MaterialManagement HeatFDTD PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src" )

# Test subproject
add_subdirectory(test)
