# Link previously defined CMake compilation "targets" together as needed

# link source targets
target_link_libraries(               interpolator2 ionize interpRateCoeff)

target_link_libraries( interp2d thrust )

target_link_libraries( flags libconfig thrust )

target_link_libraries( utils 
                       libconfig
                       thrust
                       interp2d
                       interpolator2
                       netcdf )

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
              interpolator2
              materials
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

  foreach( target IN LISTS gpu_targets  )
  # foreach( target IN LISTS gpu_targets gpu_test_targets )


    target_link_libraries( ${target} thrust CUDA::cudart )

  endforeach()

endif()
