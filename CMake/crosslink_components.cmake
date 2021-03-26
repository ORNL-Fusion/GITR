# Link previously defined CMake compilation "targets" together as needed
# link source targets
target_link_libraries( interp2d thrust )
target_link_libraries( flags libconfig thrust )
target_link_libraries( utils libconfig thrust interp2d netcdf )

# Captain! Conditionally link based on whether the GITR_USE_<comonent> clause is enabled
target_link_libraries( GITR 
                       interp2d
                       netcdf 
                       spectroscopy
                       libconfig
                       utils
                       flags
                       mpi
                       CUDA::cudart )

# link test targets
#target_link_libraries( coulomb_tests 
#                       test_utils libconfig thrust interp2d utils flags netcdf mpi )
#target_link_libraries( atomic_tests test_utils interp2d utils flags mpi )
#target_link_libraries( field_tests test_utils interp2d libconfig utils netcdf mpi fields )
#target_link_libraries( file_io_tests test_utils libconfig utils flags )
