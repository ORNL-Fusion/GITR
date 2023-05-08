cmake_minimum_required( VERSION 3.13 )

project( gitr CXX C )

set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_BUILD_TYPE "Debug" )

# point cmake to the find_package modules 
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/ )
