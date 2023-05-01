include( FindOpenMP )

if( NOT OpenMP_CXX_FOUND )

  message( FATAL_ERROR "OpenMP not found" )

endif()

include_directories( ${OpenMP_C_INCLUDE_DIRS} ${OpenMP_CXX_INCLUDE_DIRS} )
