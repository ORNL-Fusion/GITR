# Handle external dependencies

include( ExternalProject )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )

  set( suffix ".dylib" )

else()

  set( suffix ".so" )

endif()

# OpenMP
find_package( OpenMP )

if( OpenMP_CXX_FOUND )

  #set( CMAKE_CXX_FLAGS )
  message( "all the flags: ${OpenMP_CXX_FLAGS}" )
  add_compile_options( ${OpenMP_CXX_FLAGS} )

endif()

# CUDA

include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

# HDF5
include( CMake/hdf5.cmake )

# MPI

include( CMake/mpi.cmake ) # ---> creates target mpi

set( dependencies "" )

# Thrust

include( CMake/thrust.cmake ) # ---> creates target thrust

list( APPEND dependencies thrust )

# Libconfig

include( CMake/libconfig.cmake ) # ---> creates target libconfig

list( APPEND dependencies libconfig )

# NETCDF

include( CMake/netcdf.cmake ) # ---> creates target netcdf

list( APPEND dependencies netcdf )
