# Configure system for GPU support if specified

if( GITR_USE_CUDA )

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA )

  set( CUDA_SEPARABLE_COMPILATION ON )

  set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})

  # Note: CMAKE_CUDA_COMPILER must be set so that the enable_language() call succeeds
  enable_language( CUDA )

  # Later: enforce minimum requirement
  include( FindCUDAToolkit )

  include_directories( ${CUDAToolkit_INCLUDE_DIRS} )

  # Take care that this indicates the correct version
  set( CMAKE_CUDA_ARCHITECTURES 70 )

  #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -DCUDA -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA --std=c++14 -O3 --expt-relaxed-constexpr --expt-extended-lambda) #-O3 --expt-extended-lambda --expt-relaxed-constexpr -g -G --cudart shared
  #set( CMAKE_CUDA_FLAGS --expt-relaxed-constexpr )

else()

  add_compile_definitions( THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP )

endif()
