# define the order that targets must be built in to respect dependencies

# ensure that all targets in ${dependencies} built in dependencies.cmake are completely built
# before linking 
if( dependencies )

  foreach( component IN LISTS source_components test_components )

    add_dependencies( ${component} ${dependencies} )

  endforeach()

endif()

# ensure that all targets in ${source_components} are built before any test targets

foreach( component IN LISTS test_components )
  add_dependencies( ${component} ${source_components})
endforeach()

# ensure that all source targets are built before GITR
add_dependencies( GITR ${source_components} )
