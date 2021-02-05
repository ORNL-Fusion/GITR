# define singleton components, no crosslinking
set( source_components 
     efield_interp
     interp2d
     particle
     utils
     flags
     setup)

foreach( component IN LISTS source_components )

add_library( ${component} src/${component}.cpp )

target_include_directories( ${component} PRIVATE include )

endforeach()
