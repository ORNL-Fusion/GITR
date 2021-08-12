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

# CUDA
include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

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
