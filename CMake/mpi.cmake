# configure system for MPI support if specified

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

    message( FATAL_ERROR "MPI was not found - please specify MPI library location with" 
             " -DMPI_HOME=/path/to/mpi or disable MPI support with -DGITR_USE_MPI=0" )

  endif()

endif()
