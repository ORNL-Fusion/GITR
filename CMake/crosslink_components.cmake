# Link previously defined CMake compilation "targets" together as needed
# link source targets
target_link_libraries( interp2d thrust )
target_link_libraries( flags libconfig thrust )
target_link_libraries( utils libconfig thrust interp2d netcdf )
target_link_libraries( GITR interp2d netcdf libconfig utils flags mpi )

# link test targets
target_link_libraries( coulomb_tests 
                       test_utils libconfig thrust interp2d utils flags netcdf mpi )
target_link_libraries( atomic_tests test_utils interp2d utils flags mpi )
target_link_libraries( field_tests test_utils interp2d libconfig utils netcdf mpi )
target_link_libraries( file_io_tests test_utils libconfig utils flags )
