# define source components - CMake "targets" - as separate compilation components

if( USE_CUDA )
  set_source_files_properties( src/gitr.cpp PROPERTIES LANGUAGE CUDA )
endif()

add_executable( GITR src/gitr.cpp )

target_include_directories( GITR PRIVATE include )

set( source_components 
     efield_interp
     interp2d
     particle
     utils
     flags
     setup)

foreach( component IN LISTS source_components )

  add_library( ${component} src/${component}.cpp )

  if( USE_CUDA )

    set_source_files_properties( src/${component}.cpp PROPERTIES LANGUAGE CUDA )
    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )

  endif()

  target_include_directories( ${component} PUBLIC include )

endforeach()

# Add sources not in standard locations

target_sources( interp2d PUBLIC include/interp2d.hpp )
