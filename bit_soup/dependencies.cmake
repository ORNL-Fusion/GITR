include( ExternalProject )

set( prefix "${CMAKE_BINARY_DIR}/external" )

if( APPLE )

  set( suffix ".dylib" )

else()

  set( suffix ".so" )

endif()

set( CMAKE_BUILD_WITH_INSTALL_RPATH True )

# Captain! New dirty code!
include( FindHDF5 )

if( NOT HDF5_FOUND )

  message( FATAL_ERROR "Ahoy, Captain! The ship has sunk!" )

endif()

include_directories( ${HDF5_INCLUDE_DIRS} )

message( "Ahoy, Captain! ${HDF5_INCLUDE_DIRS}" )
# end dirty code

# Catch2
#include( catch.cmake )
find_package( Catch2 3 REQUIRED )

# netcdf
include( netcdf.cmake )

include( libconfig.cmake )














