# Handle external dependencies

include( ExternalProject )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )

  set( suffix ".dylib" )

else()

  set( suffix ".so" )

endif()

set( CMAKE_BUILD_WITH_INSTALL_RPATH True )

set( CMAKE_INSTALL_RPATH
     "${prefix}/libconfig_install/lib"
     "${prefix}/netcdf-c-install/lib"
     "${prefix}/netcdf-cxx4-install/lib" )

set( CMAKE_INSTALL_RPATH_USE_LINK_PATH True )

# OpenMP
find_package( OpenMP )

# Captain! Will this work with clang++ as well? brew install it and check

message( "Captain! OpenMP_C_INCLUDE_DIRS: ${OpenMP_C_INCLUDE_DIRS}" )
message( "Captain! OpenMP_CXX_INCLUDE_DIRS: ${OpenMP_CXX_INCLUDE_DIRS}" )
message( "Captain! OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
message( "Captain! OpenMP_CXX_FLAGS: ${OpenMP_C_FLAGS}")
include_directories( ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS} )

# CUDA
include( CMake/cuda.cmake ) # ---> creates target CUDA::cudart

# MPI
include( CMake/mpi.cmake ) # ---> creates target mpi

# hdf5
include( FindHDF5 )

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
