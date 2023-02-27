# Link previously defined CMake compilation "targets" together as needed

# link source targets
target_link_libraries( ionize interpRateCoeff )

target_link_libraries( interp2d thrust )

target_link_libraries( flags libconfig thrust )

target_link_libraries( fields utils )

target_link_libraries( utils 
                       libconfig
                       thrust
                       interp2d
                       netcdf )

target_link_libraries( boris thrust )

target_link_libraries( boris_data_broker boris flags )

# Captain! Does this need boris? I took it out...
target_link_libraries( cross_field_diffusion_broker
                       utils flags libconfig spectroscopy thrust
                       geometry_check config_interface )

target_link_libraries( coulomb_data_broker
                       libconfig thrust interp2d utils flags netcdf boris fields )

target_link_libraries( atomic_data_broker ionize interp2d utils flags )

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
                       netcdf 
                       spectroscopy
                       libconfig
                       utils
                       boris
                       surface_model
                       flags
                       hashGeom
                       geometry_check
                       cli11
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
endif()

# link test targets
target_link_libraries( config_interface_tests libconfig config_interface Catch2::Catch2WithMain )

target_link_libraries( coulomb_tests coulomb_data_broker
                       libconfig thrust interp2d utils flags netcdf boris fields 
                       Catch2::Catch2WithMain )

target_link_libraries( atomic_tests atomic_data_broker ionize 
                       interp2d utils flags Catch2::Catch2WithMain libconfig )

target_link_libraries( file_io_tests 
                       libconfig utils flags boris geometry_check Catch2::Catch2WithMain )

target_link_libraries( cross_field_diffusion_tests cross_field_diffusion_broker
                       utils flags libconfig boris spectroscopy thrust
                       geometry_check config_interface Catch2::Catch2WithMain )

target_link_libraries( boris_tests boris_data_broker flags libconfig utils boris slow_math Catch2::Catch2WithMain )

target_link_libraries( slow_math_tests slow_math Catch2::Catch2WithMain )

target_link_libraries( interpolator_tests Catch2::Catch2WithMain )

if( OpenMP_CXX_FOUND )

  target_link_libraries( cross_field_diffusion_tests OpenMP::OpenMP_CXX )
  target_link_libraries( coulomb_tests OpenMP::OpenMP_CXX )
  target_link_libraries( atomic_tests OpenMP::OpenMP_CXX )
  target_link_libraries( boris_tests OpenMP::OpenMP_CXX )

endif()

target_link_libraries( hdf5_tests hdf5 Catch2::Catch2WithMain )

target_link_libraries( interpolator_tests Catch2::Catch2WithMain hdf5 )

target_link_libraries( netcdf_tests Catch2::Catch2WithMain netcdf )

#add_dependencies( config_interface_tests config_interface )

#target_link_libraries( config_interface_tests libconfig config_interface Catch2::Catch2WithMain )
