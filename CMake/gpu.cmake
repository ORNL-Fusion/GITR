# create a component to encapsulate the external testing framework
# Since catch2 is a header-only library, it's target can be created as an interface
# Anything that wants to use catch2 can link to this
add_library( catch2 INTERFACE )

target_include_directories( catch2 INTERFACE 
                            test/include 
                            external/test/include )

add_library( test_utils test/src/test_utils.cpp test/include/test_utils.hpp )
target_include_directories( test_utils PUBLIC ${CMAKE_SOURCE_DIR} )

# does this propagate the includes from the catch2 interface target?
# is a target_include_directories() needed here too?
target_link_libraries( test_utils PUBLIC catch2 )

set( component 
     atomic_tests )

add_executable( ${component} test/src/${component}.cpp )
target_include_directories( ${component} PUBLIC include )

if( USE_CUDA )
set_source_files_properties( test/src/${component}.cpp PROPERTIES LANGUAGE CUDA )
set_target_properties( ${component} PROPERTIES LINKER_LANGUAGE CUDA )

# mark all files in the include directory as sources and link them
# use the dreaded glob...
file( GLOB source_files include/* test/include/* )
foreach( source_file IN LISTS source_files )
target_sources( ${component} PUBLIC ${source_file} )
set_source_files_properties( ${source_file} PROPERTIES LANGUAGE CUDA )
endforeach()
endif()
