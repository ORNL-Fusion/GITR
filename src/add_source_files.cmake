# This file adds sources in this directory to the top level cmake project
set( sources
     EfieldInterp.cpp
     Particle.cpp
     flags.cpp
     gitr.cpp
     interp2d.cpp
     setup.cpp
     utils.cpp)

# use a list operation to prepend the directory
list( TRANSFORM sources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/" OUTPUT_VARIABLE sources )

target_sources( GITR PRIVATE ${sources} )
