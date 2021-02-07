# define singleton components, no crosslinking
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

# add header files as sources where necessary - this may or may not have done anything
target_sources( interp2d PUBLIC include/interp2d.hpp )
if( USE_CUDA )
set_source_files_properties( include/interp2d.hpp PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/thermalForce.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/hashGeom.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/sortParticles.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/Surfaces.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/flags.hpp PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/interpRateCoeff.hpp PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/utils.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/recombine.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/curandInitialize.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/ionize.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/Particles.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/Particle.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/fieldLineTrace.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/ompPrint.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/h1.cuh PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/geometryCheck.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/hashGeomSheath.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/surfaceModel.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/coulombCollisions.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/testIonize.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/interpolate.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/Fields.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/history.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/testRoutine.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/testRoutineArray.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/crossFieldDiffusion.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/boris.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/testRoutineCuda.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/curandInitialize2.h PROPERTIES LANGUAGE CUDA )
set_source_files_properties( include/spectroscopy.h PROPERTIES LANGUAGE CUDA )
endif()
