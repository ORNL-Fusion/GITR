# define the order that targets must be built in to respect dependencies

# ensure that all targets in ${dependencies} built in dependencies.cmake are completely built
# before linking 
if( dependencies )

  foreach( component IN LISTS cpu_targets gpu_targets )

    add_dependencies( ${component} ${dependencies} )

  endforeach()

endif()

# ensure that test targets are built after source targets
foreach( target IN LISTS cpu_test_targets gpu_test_targets )

  add_dependencies( ${target} ${cpu_targets} ${gpu_targets} )

endforeach()

# ensure that all source targets are built before GITR
add_dependencies( GITR ${cpu_targets} ${gpu_targets} )
