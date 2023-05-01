# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

enable_testing()

set( non_gpu_test_targets
     slow_math_tests
     interpolator_tests
     config_interface_tests )

# atomic tests does not compile and is disabled
set( gpu_test_targets
     file_io_tests
     coulomb_tests 
     field_tests
     atomic_tests 
     boris_tests
     cross_field_diffusion_tests )

if( NOT GITR_USE_CUDA )

  set( non_gpu_test_targets ${non_gpu_test_targets} ${gpu_test_targets} )

endif()

foreach( component IN LISTS non_gpu_test_targets )

  add_executable( ${component} test_src/${component}.cpp )

  add_test( NAME ${component} COMMAND ${component} )

  target_include_directories( ${component} PUBLIC include test_include )

endforeach()

if( GITR_USE_CUDA )

  foreach( component IN LISTS gpu_test_targets )

    add_executable( ${component} test_src/${component}.cpp )

    add_test( NAME ${component} COMMAND ${component} )

    target_include_directories( ${component} PUBLIC include test_include )

    set_source_files_properties( test_src/${component}.cpp PROPERTIES LANGUAGE CUDA )

    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )

    set_target_properties( ${component} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON )

    target_compile_options( ${component} PUBLIC --expt-relaxed-constexpr )

  endforeach()

endif()

# Captain! Add a "ninja clean that will delete all the 
# add simple system tests
add_test( NAME system_particle_straightline 
          COMMAND GITR -c input/gitrInput.cfg
          WORKING_DIRECTORY
          "${CMAKE_BINARY_DIR}/examples/particle_trajectories/straightLine/2Dgeom" )

add_test( NAME system_particle_gyro
          COMMAND GITR -c input/gitrInput.cfg
          WORKING_DIRECTORY
          "${CMAKE_BINARY_DIR}/examples/particle_trajectories/gyroMotion" )
