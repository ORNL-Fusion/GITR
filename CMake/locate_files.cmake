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

set( dependencies "" )

# The following lines populate the "dependencies" variable

# Catch2
include( CMake/catch.cmake )

# CLI11
# include( CMake/CLI11.cmake ) # ---> creates target cli11

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

#find_package( HDF5 COMPONENTS C HL )

#if( HDF5_FOUND )
#  message( "Captain! HDF5 located!" )
#endif()

# NETCDF
include( CMake/netcdf.cmake ) # ---> creates target netcdf

# Captain! Can it operate without this? I think it probably can...
find_package( Python3 COMPONENTS Interpreter )
