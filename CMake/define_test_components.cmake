# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

enable_testing()

add_test(NAME sample_test COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a/test.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/sft_a )

set( cpu_test_targets
    )

# gpu test targets first
set( gpu_test_targets 
     boris_tests 
     cross_field_diffusion_tests
     interpolator_tests
     atomic_tests
     coulomb_tests
     file_io_tests
     )

foreach( component IN LISTS cpu_test_targets )

  add_executable( ${component} test_src/${component}.cpp )

  add_test( NAME ${component} COMMAND ${component} )

  target_include_directories( ${component} PUBLIC include test_include )

endforeach()

# compile with nvcc
if( GITR_USE_CUDA )

  # process each file - these are cpp files but link to non-cpp files
  foreach( component IN LISTS gpu_test_targets )

    add_executable( ${component} test_src/${component}.cpp )

    add_test( NAME ${component} COMMAND ${component} )

    target_include_directories( ${component} PUBLIC include test_include )

    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" ) 

    set_target_properties( ${component} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON ) 

  endforeach()


# compile gpu targets as cpu targets with host compiler
else()

  foreach( component IN LISTS gpu_test_targets )

    add_executable( ${component} test_src/${component}.cpp )

    add_test( NAME ${component} COMMAND ${component} )

    target_include_directories( ${component} PUBLIC include test_include )

  endforeach()

endif()
