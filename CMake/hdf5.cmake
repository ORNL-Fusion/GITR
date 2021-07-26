# add hdf5 if not found
include( FindHDF5 )

if( HDF5_FOUND )

  message( STATUS "HDF5 found" )

else()

  message( STATUS "HDF5 not found" )

endif()
