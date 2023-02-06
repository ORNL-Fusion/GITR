# trash

set( HDF5_PREFIX "/home/5n4/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-11.2.0/hdf5-1.10.8-tpvlximxhymutveo7jn5y2wzbhyoc2rj" )

set( NETCDF_C_PREFIX "/home/5n4/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-11.2.0/netcdf-c-4.8.1-mjek7dkgfrqg6tml7sgywdp3fjqyl6qx")

set( NETCDF_CXX_PREFIX "/home/5n4/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-11.2.0/netcdf-cxx4-4.3.1-xygkmoabzkvwyvbu5prxjx6lueywqwru" )

set( HDF5_INCLUDE_DIRS "${HDF5_PREFIX}/include" )
set( HDF5_LIBRARIES "${HDF5_PREFIX}/lib/libhdf5.so" "${HDF5_PREFIX}/lib/libhdf5_hl.so" )
set( NETCDF_INCLUDE_DIR "${NETCDF_C_PREFIX}/include" )

set( NETCDF_CXX_INCLUDE_DIR "${NETCDF_CXX_PREFIX}/include" )

set( NETCDF_CXX_LIBRARY "${NETCDF_CXX_PREFIX}/lib/libnetcdf_c++4.so" )
set( NETCDF_LIBRARY "${NETCDF_C_PREFIX}/lib/libnetcdf.so" )

# end new code 

include_directories( ${NETCDF_CXX_INCLUDE_DIR} )

include_directories( ${NETCDF_INCLUDE_DIR} )

add_library( netcdf INTERFACE )

target_include_directories( netcdf INTERFACE
                            ${NETCDF_INCLUDE_DIR}
                            ${NETCDF_CXX_INCLUDE_DIR} )

target_link_libraries( netcdf INTERFACE
                       ${NETCDF_CXX_LIBRARY} )

list( APPEND dependencies netcdf )
