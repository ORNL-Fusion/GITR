# Configure system for GPU support if specified

set( CUDA_SEPARABLE_COMPILATION ON )

set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

enable_language( CUDA )

# Captain: enforce minimum requirement
include( FindCUDAToolkit )

if( NOT CUDAToolkit_FOUND )

  message( FATAL_ERROR
  "CUDA toolkit not found: to enable, set -DCUDAToolkit_ROOT=/path/to/cuda_root."
  "The CUDA root path can be inferred from your nvcc executable: run 'which nvcc'
  to indicate where nvcc is located, take the directory portion of the path as
  the root."
  " Alternatively, to disable GPU support, set -DGITR_USE_CUDA=0 when configuring
  with CMake" )

endif()

if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )

  set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G" )

endif()

include_directories( ${CUDAToolkit_INCLUDE_DIRS} )

if( CMAKE_MINOR_VERSION VERSION_LESS 25 )
set( CMAKE_CUDA_ARCHITECTURES 80 75 70 )
else()
set( CMAKE_CUDA_ARCHITECTURES native )
endif()
