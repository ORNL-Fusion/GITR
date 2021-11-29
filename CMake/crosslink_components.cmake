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

target_link_libraries( boris thrust )

if( OpenMP_CXX_FOUND )

  target_compile_options( surface_model PUBLIC ${OpenMP_CXX_FLAGS} )
  # Captain! Is this needed? Do the compiler flags correctly handle everything? I bet they do
  target_link_libraries( surface_model PUBLIC OpenMP::OpenMP_CXX )

endif()

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
                       hashGeom
                       config_interface )

if( GITR_USE_CUDA )

  target_link_libraries( GITR 
                         CUDA::cudart )

endif()

if( GITR_USE_MPI )
  target_link_libraries( GITR mpi )
  target_link_libraries( coulomb_tests mpi )
  target_link_libraries( atomic_tests mpi )
  target_link_libraries( field_tests mpi )
endif()

# link test targets
target_link_libraries( config_interface_tests catch2 libconfig config_interface )
target_link_libraries( coulomb_tests 
                       catch2 libconfig thrust interp2d utils flags netcdf boris fields )
target_link_libraries( atomic_tests catch2 ionize interp2d utils flags )
target_link_libraries( field_tests catch2 interp2d libconfig utils netcdf fields boris )
target_link_libraries( file_io_tests catch2 libconfig utils flags boris )
target_link_libraries( cross_field_diffusion_tests 
                       catch2 utils flags libconfig boris spectroscopy )
target_link_libraries( boris_tests catch2 flags libconfig utils boris )
target_link_libraries( interpolator_tests catch2 )
