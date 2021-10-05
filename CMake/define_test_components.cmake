# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

enable_testing()

set( cpu_test_targets
     config_interface_tests
     interpolator_tests )

# atomic tests does not compile and is disabled
set( gpu_test_targets
     file_io_tests
     coulomb_tests 
     field_tests
     atomic_tests 
     boris_tests
     cross_field_diffusion_tests )

if( NOT GITR_USE_CUDA )

  set( cpu_test_targets ${cpu_test_targets} ${gpu_test_targets} )

endif()

foreach( component IN LISTS cpu_test_targets )

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

    set_target_properties( ${component} PROPERTIES
                           COMPILE_FLAGS "-dc"
                           CUDA_RESOLVE_DEVICE_SYMBOLS ON )

    target_compile_options( ${component} PUBLIC --expt-relaxed-constexpr )

  endforeach()

endif()
