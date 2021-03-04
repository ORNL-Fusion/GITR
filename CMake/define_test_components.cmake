# define test components - CMake "targets" - as separate compilation components

# create a component to encapsulate the external testing framework
# Since catch2 is a header-only library, it's target can be created as an interface target
add_library( catch2 INTERFACE )

target_include_directories( catch2 INTERFACE 
                            test_include )

add_library( test_utils test_src/test_utils.cpp test_include/test_utils.hpp )
target_include_directories( test_utils PUBLIC ${CMAKE_SOURCE_DIR} )

target_link_libraries( test_utils PUBLIC catch2 )

set( test_components 
     atomic_tests
     coulomb_tests
     field_tests
     file_io_tests )

foreach( component IN LISTS test_components )

  add_executable( ${component} test_src/${component}.cpp )

  if( USE_CUDA )

    set_source_files_properties( test_src/${component}.cpp PROPERTIES LANGUAGE CUDA )

    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )

  endif()

  target_include_directories( ${component} PUBLIC include test_include )

endforeach()
