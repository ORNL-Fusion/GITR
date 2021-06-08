# copies over the file and sets destination_path
macro( generate_testing_file test_file_name )

set( source_path "${CMAKE_SOURCE_DIR}/${test_file_name}" )
set( destination_path "${CMAKE_BINARY_DIR}/${test_file_name}" )
configure_file( ${source_path} ${destination_path} )

endmacro()

generate_testing_file( "test_data/test_config.cfg" )
set( CONFIG_INTERFACE_UNIT_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/test_coulomb.cfg" )
set( COULOMB_UNIT_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/ionize.cfg" )
set( FIELD_UNIT_TEST_FILE_0 ${destination_path} )

generate_testing_file( "test_data/ionize1.cfg" )
set( FIELD_UNIT_TEST_FILE_1 ${destination_path} )

# Configure the header file that will contain these strings

set( source_name "test_include/test_data_filepath.hpp.in" )

set( source_path "${CMAKE_SOURCE_DIR}/${source_name}" )

set( destination_name "test_include/test_data_filepath.hpp" )

set( destination_path "${CMAKE_SOURCE_DIR}/${destination_name}" )

configure_file( ${source_path} ${destination_path} )
