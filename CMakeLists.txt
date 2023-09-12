cmake_minimum_required( VERSION 3.13 )

project( gitr CXX C )

set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_BUILD_TYPE "Release" )

# point cmake to the find_package modules 
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/ )

# CMake options configurable with -D<option_string>
include( CMake/user_options.cmake )

# Handle external dependencies
include( CMake/dependencies.cmake )

# define singleton source components - change to "targets"
include( CMake/define_source_components.cmake )

## AD: comment this for now
# # define singleton test components
# include( CMake/define_test_components.cmake )

# configure testing data
include( CMake/configure_test_data.cmake )

# Make sure everything is built before attempting to link
include( CMake/enforce_build_order.cmake )

# Link targets together for final executable
include( CMake/crosslink_components.cmake )
