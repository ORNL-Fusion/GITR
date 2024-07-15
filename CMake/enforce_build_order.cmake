# define the order that targets must be built in to respect dependencies

# ensure that all targets in ${dependencies} built in dependencies.cmake are completely built
# before linking 
if( dependencies )

  foreach( component IN LISTS non_gpu_targets gpu_targets )

    add_dependencies( ${component} ${dependencies} )

  endforeach()

endif()

add_dependencies( boris particle_diagnostics )
# ensure that all source targets are built before GITR
add_dependencies( GITR ${non_gpu_targets} ${gpu_targets} )

add_dependencies( boris_data_broker boris flags )

if( GITR_USE_CUDA )

  add_dependencies( geometry_check atomic_add_1 )
  add_dependencies( spectroscopy atomic_add_1 )

endif()
