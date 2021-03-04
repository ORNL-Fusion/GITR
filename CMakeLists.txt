cmake_minimum_required(VERSION 3.13)

project(gitr)

set( CMAKE_CXX_STANDARD 17 )

# point cmake to the find_package modules 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)

# CMake options configurable with -D<option_string>
include( CMake/user_options.cmake )

# Configures system according to user options
include( CMake/system_options.cmake )

# Handle external dependencies
include( CMake/dependencies.cmake )

# CPU mode
if( NOT GITR_USE_CUDA )

  # define singleton source components - change to "targets"
  include( CMake/define_source_components.cmake )

  # define singleton test components
  include( CMake/define_test_components.cmake )

  # Make sure everything is built before attempting to link
  include( CMake/enforce_build_order.cmake )

  include( CMake/crosslink_components.cmake )

# GPU mode not yet enabled
elseif( GITR_USE_CUDA )

  message( STATUS "GPU mode not enabled" )

endif()
