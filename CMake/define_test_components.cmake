# define test components - CMake "targets" - as separate compilation components

# create a component to encapsulate the external testing framework
# Since catch2 is a header-only library, it's target can be created as an interface target
add_library( catch2 INTERFACE )

target_include_directories( catch2 INTERFACE 
                            test_include )

add_library( test_utils test_src/test_utils.cpp test_include/test_utils.hpp )

target_include_directories( test_utils PUBLIC ${CMAKE_SOURCE_DIR} )

target_link_libraries( test_utils PUBLIC catch2 )

set( cpu_test_targets
     file_io_tests )

set( gpu_test_targets
     coulomb_tests 
     field_tests
     atomic_tests )

if( NOT GITR_USE_CUDA )

  set( cpu_test_targets ${cpu_test_targets} ${gpu_test_targets} )

endif()

foreach( component IN LISTS cpu_test_targets )

  add_executable( ${component} test_src/${component}.cpp )

  target_include_directories( ${component} PUBLIC include test_include )

endforeach()

if( GITR_USE_CUDA )

  foreach( component IN LISTS gpu_test_targets )

    add_executable( ${component} test_src/${component}.cpp )
    target_include_directories( ${component} PUBLIC include test_include )
    set_source_files_properties( test_src/${component}.cpp PROPERTIES LANGUAGE CUDA )
    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )
    target_compile_options( ${component} PUBLIC --expt-relaxed-constexpr )

  endforeach()

endif()
