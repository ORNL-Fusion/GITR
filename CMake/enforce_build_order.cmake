# define the order that targets must be built in to respect dependencies

# ensure that all targets in ${dependencies} built in dependencies.cmake are completely built
# before linking 
if( dependencies )

  foreach( component IN LISTS non_gpu_targets gpu_targets )

    add_dependencies( ${component} ${dependencies} )

  endforeach()

endif()

foreach( target IN LISTS cpu_test_targets gpu_test_targets )

  add_dependencies( ${target} ${cpu_targets} ${gpu_targets} ${dependencies} ${gpu_broker_targets} )

endforeach()

# ensure that all source targets are built before GITR
add_dependencies( GITR ${non_gpu_targets} ${gpu_targets} )

add_dependencies( boris_data_broker boris flags )

add_dependencies( cross_field_diffusion_broker ${gpu_targets} ${non_gpu_targets} ${dependencies} )

add_dependencies( atomic_data_broker ${gpu_targets} ${non_gpu_targets} ${dependencies} )

add_dependencies( coulomb_data_broker ${gpu_targets} ${non_gpu_targets} ${dependencies} )

add_dependencies( fields utils )
