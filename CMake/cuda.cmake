# Configure system for GPU support if specified

if( GITR_USE_CUDA )

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA )

  set( CUDA_SEPARABLE_COMPILATION ON )

  set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

  # Note: CMAKE_CUDA_COMPILER must be set so that the enable_language() call succeeds
  enable_language( CUDA )

  # Later: enforce minimum requirement

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

  include_directories( ${CUDAToolkit_INCLUDE_DIRS} )

  set( CMAKE_CUDA_ARCHITECTURES 70 )

else()

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP )

endif()
