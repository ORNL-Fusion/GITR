# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

enable_testing()

set( non_gpu_test_targets
     slow_math_tests
     config_interface_tests )

# Captain! Next, you simply have to go through and re-activate the rest of the tests.
# once all the tests are reactivated, test them all out in CPU mode too. Then add in
# your new tests. No documentation, do that tomorrow
# atomic tests does not compile and is disabled
set( gpu_test_targets
     file_io_tests
     coulomb_tests 
     atomic_tests 
     cross_field_diffusion_tests
     boris_tests )

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

    #set_source_files_properties( test_src/${component}.cpp PROPERTIES LANGUAGE CXX )
    #set_source_files_properties( test_src/${component}.cpp PROPERTIES LINKER_LANGUAGE CUDA )

    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )

    set_target_properties( ${component} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON )

    #target_compile_options( ${component} PUBLIC --expt-relaxed-constexpr )

  endforeach()

endif()

add_test( NAME system_particle_straightline 
          COMMAND GITR -c input/gitrInput.cfg
          WORKING_DIRECTORY
          "${CMAKE_BINARY_DIR}/examples/particle_trajectories/straightLine/2Dgeom" )

add_test( NAME system_particle_gyro
          COMMAND GITR -c input/gitrInput.cfg
          WORKING_DIRECTORY
          "${CMAKE_BINARY_DIR}/examples/particle_trajectories/gyroMotion" )

# Captain! New code
add_executable( hdf5_tests test_src/hdf5_tests.cpp )

add_test( NAME hdf5_tests COMMAND hdf5_tests )

target_include_directories( hdf5_tests PUBLIC include test_include )

add_executable( interpolator_tests test_src/interpolator_tests.cpp )

add_test( NAME interpolator_tests COMMAND interpolator_tests )

target_include_directories( interpolator_tests PUBLIC include test_include )

add_executable( netcdf_tests test_src/netcdf_tests.cpp )

add_test( NAME netcdf_tests COMMAND netcdf_tests )

target_include_directories( netcdf_tests PUBLIC include test_include )
