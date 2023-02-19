# Configure system for GPU support if specified

set( CUDA_SEPARABLE_COMPILATION ON )

set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

enable_language( CUDA )

set( CUDAToolkit_ROOT "/usr/local/cuda/include" CACHE STRING "" FORCE )

# Captain: enforce minimum requirement

# Captaaaiiiin! Looks like you may have to install the cuda toolkit on your system...
# Orrrrr can you do it in the container?? Perhaps... we shall see
# Captain! What if it's just apt install cuda... try that one out now...
include( FindCUDAToolkit )

if( CUDAToolkit_FOUND )

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

set( CMAKE_CUDA_ARCHITECTURES 70 )
