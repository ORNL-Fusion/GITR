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
     "${prefix}/netcdf-cxx4-install/lib"
     "${prefix}/CMake-hdf5-1.13.0/build/bin" )

# OpenMP
find_package( OpenMP )

include_directories( ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS} )

# CLI11
include( CMake/CLI11.cmake )

list( APPEND dependencies cli11 )

# hdf5
# Captain! Get rid of the "find" modules - those are really outdated
include( CMake/hdf5.cmake )

list( APPEND dependencies hdf5 )

# CUDA
include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

# MPI
include( CMake/mpi.cmake ) # ---> creates target mpi

# Thrust
include( CMake/thrust.cmake ) # ---> creates target thrust

list( APPEND dependencies thrust )

# Libconfig
include( CMake/libconfig.cmake ) # ---> creates target libconfig

list( APPEND dependencies libconfig )

# NETCDF
include( CMake/netcdf.cmake ) # ---> creates target netcdf

list( APPEND dependencies netcdf )
