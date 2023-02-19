include( CTest )

enable_testing()

add_executable( hdf5_tests hdf5_tests.cpp )

add_test( NAME hdf5_tests COMMAND hdf5_tests )

add_executable( interpolator_tests interpolator_tests.cpp )

add_test( NAME interpolator_tests COMMAND interpolator_tests )

add_executable( netcdf_tests netcdf_tests.cpp )

add_test( NAME netcdf_tests COMMAND netcdf_tests )

add_executable( config_interface_tests config_interface_tests.cpp )

add_test( NAME config_interface_tests COMMAND config_interface_tests )

add_library( config_interface config_interface.cpp )
