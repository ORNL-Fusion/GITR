# required first two lines of CMake file
cmake_minimum_required( VERSION 3.13 )

project( gitr CXX C )

set( CMAKE_INSTALL_RPATH_USE_LINK_PATH True )

# CMake options configurable with -D<option_string>
include( CMake/compile_time_options.cmake )

# Handle external dependencies
include( CMake/dependencies.cmake )

# define singleton source components - change to "targets"
include( CMake/define_source_components.cmake )

# define testing related components
include( CMake/define_test_components.cmake )

# configure testing data
include( CMake/configure_test_data.cmake )

# Make sure everything is built before attempting to link
include( CMake/enforce_build_order.cmake )

# Link targets together for final executable
include( CMake/crosslink_components.cmake )
