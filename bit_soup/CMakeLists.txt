# Project for hdf5 file I/O tests with multidimensional data arrays
cmake_minimum_required( VERSION 3.22 )

set( CMAKE_VERBOSE_MAKEFILE TRUE )

project( hdf5_file_io C CXX )

set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_BUILD_TYPE "Release" )

include( dependencies.cmake )

include( define_test_components.cmake )

include( enforce_build_order.cmake )

include( crosslink_components.cmake )
