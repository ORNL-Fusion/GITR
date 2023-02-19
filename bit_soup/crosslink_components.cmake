
target_link_libraries( hdf5_tests Catch2::Catch2WithMain  hdf5::hdf5 )

target_link_libraries( interpolator_tests Catch2::Catch2WithMain  hdf5::hdf5 )

target_link_libraries( netcdf_tests Catch2::Catch2WithMain netcdf )

add_dependencies( config_interface_tests config_interface )

target_link_libraries( config_interface_tests libconfig config_interface Catch2::Catch2WithMain )
