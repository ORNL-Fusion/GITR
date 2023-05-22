# Link previously defined CMake compilation "targets" together as needed

target_link_libraries( efield_interp libconfig_cxx libconfig_c )

# link source targets
target_link_libraries( ionize interpRateCoeff )

target_link_libraries( interp2d thrust )

target_link_libraries( flags libconfig_cxx libconfig_c thrust )

target_link_libraries( utils 
                       libconfig_cxx
                       libconfig_c
                       thrust
                       interp2d
                       netcdf_cxx
                       netcdf_c )

target_link_libraries( boris thrust )

if( OpenMP_CXX_FOUND )

  target_compile_options( surface_model PUBLIC ${OpenMP_CXX_FLAGS} )
  target_link_libraries( surface_model OpenMP::OpenMP_CXX )

  target_compile_options( spectroscopy PUBLIC ${OpenMP_CXX_FLAGS} )
  target_link_libraries( spectroscopy OpenMP::OpenMP_CXX )

  target_compile_options( geometry_check PUBLIC ${OpenMP_CXX_FLAGS} )
  target_link_libraries( geometry_check OpenMP::OpenMP_CXX )

endif()

target_link_libraries( geometry_check boris )

# Improvement: Conditionally link based on whether the GITR_USE_<component> clause is enabled
target_link_libraries( GITR 
                       ionize
                       interp2d
                       ${netcdf_cxx_shared_lib}
                       ${netcdf_c_shared_lib}
                       spectroscopy
                       libconfig_cxx
                       libconfig_c
                       utils
                       boris
                       surface_model
                       flags
                       hashGeom
                       geometry_check
                       config_interface )

if( GITR_USE_CUDA )

  target_link_libraries( GITR 
                         CUDA::cudart )

  foreach( target IN LISTS gpu_targets gpu_test_targets )

    target_link_libraries( ${target} thrust CUDA::cudart )

  endforeach()

endif()

if( GITR_USE_MPI )
  target_link_libraries( GITR mpi )
  target_link_libraries( coulomb_tests mpi )
  target_link_libraries( atomic_tests mpi )
  target_link_libraries( field_tests mpi )
endif()
