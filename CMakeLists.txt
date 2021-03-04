cmake_minimum_required(VERSION 3.13)
project(gitr)

set( CMAKE_CXX_STANDARD 17 )

# point cmake to the find_package modules 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMake/)

include( CMake/define_options.cmake )

include( CMake/dependencies.cmake )

# code path for CPU mode
if( NOT GITR_USE_CUDA )

# for now, hardcode compile definitions CPU build only.
# move this into system_config_backend
include( CMake/configure_hardware.cmake )

# define singleton source components - change to "targets"
include( CMake/define_source_components.cmake )

# define singleton test components
include( CMake/define_test_components.cmake )

# Make sure everything is built before attempting to link
include( CMake/enforce_build_order.cmake )

# not yet implemented link the targets together as needed
include( CMake/crosslink_components.cmake )

# not yet implemented - this will call "add_test"
include( CMake/register_tests.cmake )

# code path for GPU mode
elseif( GITR_USE_CUDA )
message( STATUS "GPU mode not enabled" )
endif()
