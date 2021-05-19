# Link previously defined CMake compilation "targets" together as needed
# link source targets
target_link_libraries( ionize interpRateCoeff )

target_link_libraries( interp2d thrust )

target_link_libraries( flags libconfig thrust )

target_link_libraries( utils 
                       libconfig
                       thrust
                       interp2d
                       netcdf )

# Improvement: Conditionally link based on whether the GITR_USE_<component> clause is enabled
target_link_libraries( GITR 
                       ionize
                       interp2d
                       netcdf 
                       spectroscopy
                       libconfig
                       utils
                       boris
                       surface_model
                       flags
                       hashGeom )

if( GITR_USE_CUDA )

  target_link_libraries( GITR 
                         CUDA::cudart )

endif()

if( GITR_USE_MPI )
  target_link_libraries( GITR mpi )
endif()

# link test targets
#target_link_libraries( coulomb_tests 
#                       test_utils libconfig thrust interp2d utils flags netcdf mpi )
#target_link_libraries( atomic_tests test_utils interp2d utils flags mpi )
#target_link_libraries( field_tests test_utils interp2d libconfig utils netcdf mpi fields )
#target_link_libraries( file_io_tests test_utils libconfig utils flags )
