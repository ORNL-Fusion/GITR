# CPU mode only - ensure that all targets in ${dependencies} are built before linking 
if( dependencies )

  foreach( component IN LISTS source_components test_components )
  # does this really have to be dereferenced?
  add_dependencies( ${component} ${dependencies} )
  endforeach()

endif()

# ensure that all targets in ${source_components} are built before any test targets
foreach( component IN LISTS test_components )
  add_dependencies( ${component} ${source_components})
endforeach()

# ensure that all source targets are built before GITR
add_dependencies( GITR ${source_components} )
