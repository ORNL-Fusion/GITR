# copies over the file and sets destination_path
macro( generate_testing_file test_file_name )

  set( source_path "${CMAKE_SOURCE_DIR}/${test_file_name}" )

  set( destination_path "${CMAKE_BINARY_DIR}/${test_file_name}" )

  configure_file( ${source_path} ${destination_path} COPYONLY )

endmacro()

generate_testing_file( "test_data/Efield.txt" )
set( E_FIELD_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/getE.cfg" )
set( GET_E_TEST_FILE  ${destination_path} )

generate_testing_file( "test_data/surface_model.cfg" )
set( SURFACE_MODEL_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/test_config.cfg" )
set( CONFIG_INTERFACE_UNIT_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/test_coulomb.cfg" )
set( COULOMB_UNIT_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/ionize.cfg" )
set( FIELD_UNIT_TEST_FILE_0 ${destination_path} )

generate_testing_file( "test_data/ionize1.cfg" )
set( FIELD_UNIT_TEST_FILE_1 ${destination_path} )

generate_testing_file( "test_data/cross_field_geometry.cfg" )
set( CROSS_FIELD_GEOM_FILE ${destination_path} )

generate_testing_file( "test_data/cross_field_geometry2.cfg" )
set( CROSS_FIELD_GEOM_FILE_1 ${destination_path} )

generate_testing_file( "test_data/boris_config.cfg" )
set( BORIS_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/file.cfg" )
set( DATA_TYPES_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/geom_test.cfg" )
set( GEOM_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/flat_line.cfg" )
set( FLAT_LINE_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/positive_slope.cfg" )
set( POSITIVE_SLOPE_TEST_FILE ${destination_path} )

generate_testing_file( "test_data/ADAS_Rates_W.nc" )
set( ADAS_TEST_FILE ${destination_path} )

# Configure the header file that will contain these strings

set( source_name "test_include/test_data_filepath.hpp.in" )

set( source_path "${CMAKE_SOURCE_DIR}/${source_name}" )

set( destination_name "test_include/test_data_filepath.hpp" )

set( destination_path "${CMAKE_SOURCE_DIR}/${destination_name}" )

configure_file( ${source_path} ${destination_path} )

# This .nc file wil be removed with the netcdf dependency in the future
set( test_file_name "test_data/netcdf_file_py.nc" )
set( source_path "${CMAKE_SOURCE_DIR}/${test_file_name}" )
set( destination_path "${CMAKE_BINARY_DIR}/${test_file_name}" )
configure_file( ${source_path} ${destination_path} COPYONLY )

# copy examples folder to the binary directory
set( source_path "${CMAKE_SOURCE_DIR}/examples/" )
set( destination_path "${CMAKE_BINARY_DIR}/examples" )
file( COPY ${source_path} DESTINATION ${destination_path} )
