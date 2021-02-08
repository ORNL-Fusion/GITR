# define singleton components, no crosslinking
add_executable( GITR src/gitr.cpp )
if( USE_CUDA )
set_source_files_properties( src/gitr.cpp PROPERTIES LANGUAGE CUDA )
endif()

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

target_include_directories( ${component} PRIVATE include )

endforeach()

# add header files as sources where necessary - this may or may not have done anything
target_sources( interp2d PUBLIC include/interp2d.hpp )
if( USE_CUDA )
set_source_files_properties( include/interp2d.cpp PROPERTIES LANGUAGE CUDA )
endif()
