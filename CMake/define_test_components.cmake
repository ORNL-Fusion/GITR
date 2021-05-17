# define test components - CMake "targets" - as separate compilation components
include( CTest ) 

# create a component to encapsulate the external testing framework
# Since catch2 is a header-only library, it's target can be created as an interface target
add_library( catch2 INTERFACE )

target_include_directories( catch2 INTERFACE 
                            test_include )

add_library( test_utils test_src/test_utils.cpp test_include/test_utils.hpp )

target_include_directories( test_utils PUBLIC ${CMAKE_SOURCE_DIR} )

target_link_libraries( test_utils PUBLIC catch2 )

# Captain! Disable tests here 
set( cpu_test_targets
     config_interface_tests
     # file_io_tests
     )

set( gpu_test_targets
     # coulomb_tests 
     # field_tests
     # atomic_tests
     )

if( NOT GITR_USE_CUDA )

  set( cpu_test_targets ${cpu_test_targets} ${gpu_test_targets} )

endif()

foreach( component IN LISTS cpu_test_targets )

  add_executable( ${component} test_src/${component}.cpp )

  add_test( NAME config_interface_tests COMMAND config_interface_tests )

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

# copy testing data to the build directory
set( libconfig_unit_test_file_name "test_data/test_config.cfg" )

set( libconfig_unit_test_file_source_path
     "${CMAKE_SOURCE_DIR}/${libconfig_unit_test_file_name}" )

set( libconfig_unit_test_file_destination_path 
     "${CMAKE_BINARY_DIR}/${libconfig_unit_test_file_name}" )

# For the configured file below
set( LIBCONFIG_UNIT_TEST_FILE ${libconfig_unit_test_file_destination_path} )

set( test_data_filepath_header_source_name "test_include/test_data_filepath.hpp.in" )

set( test_data_filepath_header_source_path
     "${CMAKE_SOURCE_DIR}/${test_data_filepath_header_source_name}" )

set( test_data_filepath_header_destination_name "test_include/test_data_filepath.hpp" )

set( test_data_filepath_header_destination_path
     "${CMAKE_SOURCE_DIR}/${test_data_filepath_header_destination_name}" )

configure_file( ${test_data_filepath_header_source_path} 
                ${test_data_filepath_header_destination_path} )

configure_file( ${libconfig_unit_test_file_source_path}
                ${libconfig_unit_test_file_destination_path} )
