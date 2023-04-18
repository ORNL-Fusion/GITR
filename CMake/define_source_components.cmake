# define source components - CMake "targets" - as separate compilation components

# serial modules

add_executable( GITR src/gitr.cpp )

# extra compile options for CUDA
if( GITR_USE_CUDA )

  set_source_files_properties( src/gitr.cpp PROPERTIES LANGUAGE CUDA )

  set_target_properties( GITR PROPERTIES LINKER_LANGUAGE CUDA )

  set_property(TARGET GITR PROPERTY CUDA_SEPARABLE_COMPILATION ON)

  target_compile_options( GITR PRIVATE --expt-relaxed-constexpr )

endif()

target_include_directories( GITR PUBLIC include )

# CPU-only targets
set( non_gpu_targets
     efield_interp
     particle
     utils
     flags
     config_interface
     slow_math
     interpolator
     setup)

# conditionally compile as GPU targets
set( gpu_targets
     surface_model
     interp2d
     interpRateCoeff
     ionize
     boris
     hashGeom
     fields
     spectroscopy
     geometry_check )

set( gpu_broker_targets
     boris_data_broker
     atomic_data_broker
     coulomb_data_broker
     cross_field_diffusion_broker )

# Captain! New code:
# compile all cpu targets first
foreach( component IN LISTS non_gpu_targets )

  add_library( ${component} src/${component}.cpp )

  target_include_directories( ${component} PUBLIC include test_include )

endforeach()

if( GITR_USE_CUDA )

  foreach( component IN LISTS gpu_targets gpu_broker_targets )

    add_library( ${component} src/${component}.cpp )

    set_source_files_properties( src/${component}.cpp PROPERTIES LANGUAGE CUDA )

    set_target_properties( ${component} PROPERTIES COMPILE_FLAGS "-dc" )

    target_include_directories( ${component} PUBLIC include test_include )

    target_compile_options( ${component} PRIVATE --expt-relaxed-constexpr )

  endforeach()

else()

  foreach( component IN LISTS gpu_targets gpu_broker_targets )

    add_library( ${component} src/${component}.cpp )

    target_include_directories( ${component} PUBLIC include test_include )

  endforeach()

endif()
