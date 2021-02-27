# Handle external dependencies

include( ExternalProject )

set( dependencies "" )

set( prefix "${CMAKE_CURRENT_SOURCE_DIR}/external" )

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

  else()

    message( FATAL_ERROR "MPI was not found" )

  endif()

endif()
