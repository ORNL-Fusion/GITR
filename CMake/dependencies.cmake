# Handle external dependencies
include( ExternalProject )

set( dependencies "" )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )

  set( suffix ".dylib" )

else()

  set( suffix ".so" )

endif()

set( CMAKE_BUILD_WITH_INSTALL_RPATH True )

set( CMAKE_INSTALL_RPATH
     "${CMAKE_BINARY_DIR}"
     "${prefix}/libconfig_install/lib"
     "${prefix}/netcdf-c-install/lib"
     "${prefix}/netcdf-cxx4-install/lib" )

# OpenMP
include( FindOpenMP )

include_directories( ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS} )

# CLI11
include( CMake/CLI11.cmake ) # ---> creates target cli11

list( APPEND dependencies cli11 )

# hdf5
#include( FindHDF5 )
find_package( HDF5 COMPONENTS C HL )

message( "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}" )
message( "HDF5_INCLUDE_DIRS: ${HDF5_LIBRARIES}" )

# CUDA
include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

# MPI
include( CMake/mpi.cmake ) # ---> creates target mpi

# Thrust
include( CMake/thrust.cmake ) # ---> creates target thrust

#list( APPEND dependencies thrust )

# Libconfig
include( CMake/libconfig.cmake ) # ---> creates target libconfig

list( APPEND dependencies libconfig )

# NETCDF
include( CMake/netcdf.cmake ) # ---> creates target netcdf

list( APPEND dependencies netcdf )
