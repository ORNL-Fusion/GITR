# This file adds sources in this directory to the top level cmake project
set( sources
     EfieldInterp.cpp
     Particle.cpp
     flags.cpp
     gitr.cpp
     interp2d.cpp
     setup.cpp
     utils.cpp)

target_sources( GITR PRIVATE ${sources} )
