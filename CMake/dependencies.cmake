# Handle external dependencies

include( ExternalProject )

set( dependencies "" )

set( prefix "${CMAKE_BINARY_DIR}/external" )

# Thrust

include( CMake/thrust.cmake ) # --> creates target thrust

list( APPEND dependencies thrust )

# Libconfig

include( CMake/libconfig.cmake ) # --> creates target libconfig

list( APPEND dependencies libconfig )

# NETCDF

include( CMake/netcdf.cmake ) # --> creates target netcdf

list( APPEND dependencies netcdf )

# MPI

if( GITR_USE_MPI )
  include( FindMPI )

  if(MPI_FOUND)

    add_library( mpi INTERFACE )
    target_include_directories( mpi INTERFACE ${MPI_CXX_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS})
    target_link_libraries( mpi INTERFACE ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} )

    set( MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_DIRS} CACHE PATH "" FORCE )
    set( MPI_C_INCLUDE_DIRS ${MPI_C_INCLUDE_DIRS} CACHE PATH "" FORCE )
    set( MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} CACHE FILEPATH "" FORCE )
    set( MPI_C_LIBRARIES ${MPI_C_LIBRARIES} CACHE FILEPATH "" FORCE )

  else()

    message( FATAL_ERROR "MPI was not found" )

  endif()

endif()
