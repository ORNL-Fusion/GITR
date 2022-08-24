# Only OpenMP or Cuda can be specified
if( GITR_USE_CUDA AND GITR_USE_OPENMP )

  message( FATAL_ERROR "Both GITR_USE_CUDA and GITR_USE_OPENMP are set. Please select one" )

endif()

# Handle external dependencies
include( ExternalProject )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )

  set( suffix ".dylib" )

else()

  set( suffix ".so" )

endif()

set( CMAKE_BUILD_WITH_INSTALL_RPATH True )

# ensure shared dependency libs are discoverable at load-time
set( CMAKE_INSTALL_RPATH
     "${CMAKE_BINARY_DIR}"
     "${prefix}/libconfig_install/lib"
     "${prefix}/netcdf-c-install/lib"
     "${prefix}/netcdf-cxx4-install/lib" )

# Captain! Rename to "common dependency set" or something more descriptive
set( dependencies "" )

# The following lines populate the "dependencies" variable

# CLI11
include( CMake/CLI11.cmake ) # ---> creates target cli11

# HDF5
find_package( HDF5 COMPONENTS C HL )

# CUDA
if( GITR_USE_CUDA )

  include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

# OpenMP
elseif( GITR_USE_OPENMP )

  include( CMake/openmp.cmake ) # ---> creates target OpenMP

endif()

# Thrust
include( CMake/thrust.cmake ) # ---> creates target thrust

# MPI 
if( GITR_USE_MPI )

  include( CMake/mpi.cmake ) # ---> creates target mpi

endif()

# Libconfig
include( CMake/libconfig.cmake ) # ---> creates target libconfig

# NETCDF
include( CMake/netcdf.cmake ) # ---> creates target netcdf

# Catch2
include( CMake/catch.cmake )
