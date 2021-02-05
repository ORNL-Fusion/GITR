# create a component to encapsulate the external testing framework
# Since catch2 is a header-only library, it's target can be created as an interface
# Anything that wants to use catch2 can link to this
add_library( catch2 INTERFACE )

target_include_directories( catch2 INTERFACE 
                            test/include 
                            external/test/include )

add_library( test_utils test/src/test_utils.cpp )

# does this propagate the includes from the catch2 interface target?
# is a target_include_directories() needed here too?
target_link_libraries( test_utils PUBLIC catch2 )

set( test_components 
     atomic_tests
     coulomb_tests
     field_tests
     file_io_tests )

foreach( component IN LISTS test_components )

add_library( ${component} test/src/${component}.cpp )

target_include_directories( ${component} PRIVATE include )

endforeach()

# Each component will get the include directories from catch2 when they link against it
# since it is an interface
