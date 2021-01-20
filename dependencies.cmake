# download dependencies and set targets
include( ExternalProject )

# put dependency targets here
set( dependencies "" )

find_package(Thrust)
if( NOT THRUST_FOUND )

ExternalProject_Add( thrust
                     PREFIX "/thrust" 
                     GIT_REPOSITORY "git@github.com:NVIDIA/thrust.git" 
                     CONFIGURE_COMMAND ""
                     BUILD_COMMAND ""
                     INSTALL_COMMAND "" )

set( THRUST_INCLUDE_DIR "/thrust/src/thrust" )
find_package(Thrust REQUIRED)
list( APPEND dependencies thrust )
endif()

find_package(LibConfig)
if( NOT LIBCONFIG_FOUND )

ExternalProject_Add( libconfig
                     PREFIX "/libconfig"
                     GIT_REPOSITORY "git@github.com:hyperrealm/libconfig.git" 
                     INSTALL_COMMAND "" )
set( LIBCONFIG_INCLUDE_DIR "/libconfig/src/libconfig/lib" )
set( LIBCONFIG_LIBRARY "/libconfig/src/libconfig-build/out/libconfig.so" )
set( LIBCONFIGPP_INCLUDE_DIR "/libconfig/src/libconfig/lib" )
set( LIBCONFIGPP_LIBRARY "/libconfig/src/libconfig-build/out/libconfig++.so" )
find_package(LibConfig REQUIRED)
list( APPEND dependencies libconfig )
endif()

find_package(NetCDF COMPONENTS CXX)
if( NOT NETCDF_FOUND )
ExternalProject_Add( netcdf-c
                     PREFIX "/netcdf-c"
                     GIT_REPOSITORY "git@github.com:Unidata/netcdf-c.git"
                     CMAKE_ARGS "-DENABLE_DAP=OFF" )

ExternalProject_Add( netcdf-cxx4
                     PREFIX "/netcdf-cxx4"
                     GIT_REPOSITORY "git@github.com:Unidata/netcdf-cxx4.git" )

add_dependencies( netcdf-cxx4 netcdf-c )
set( NETCDF_LIBRARY "/usr/local/lib/libnetcdf.so" )
set( NETCDF_CXX_INCLUDE_DIR "/netcdf-cxx4/src/netcdf-cxx4/cxx4" )
set( NETCDF_CXX_LIBRARY "/usr/local/lib/libnetcdf-cxx4.so" )
set( NETCDF_INCLUDE_DIR "/netcdf-c/src/netcdf-c/include" )
find_package(NetCDF COMPONENTS CXX REQUIRED)
list( APPEND dependencies netcdf-c netcdf-cxx4 )
endif()
