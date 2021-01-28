# download dependencies and set targets
include( ExternalProject )

# put dependency targets here
set( dependencies "" )

find_package(Thrust)
if( NOT THRUST_FOUND )

set( prefix "${CMAKE_CURRENT_SOURCE_DIR}/external" )

message( "Ahoy, Captain! Downloading thrust..." )
ExternalProject_Add( thrust
                     PREFIX "${prefix}/thrust" 
                     GIT_REPOSITORY "git@github.com:NVIDIA/thrust.git" 
                     CONFIGURE_COMMAND ""
                     BUILD_COMMAND ""
                     INSTALL_COMMAND "" )

set( THRUST_INCLUDE_DIR "${prefix}/thrust/src/thrust" )
find_package(Thrust REQUIRED)
list( APPEND dependencies thrust )
endif()


find_package(LibConfig)
if( NOT LIBCONFIG_FOUND )

message( "Ahoy, Captain! Downloading libconfig..." )
ExternalProject_Add( libconfig
                     PREFIX "${prefix}/libconfig"
                     GIT_REPOSITORY "git@github.com:hyperrealm/libconfig.git" 
                     INSTALL_COMMAND "" )
set( LIBCONFIG_INCLUDE_DIR "${prefix}/libconfig/src/libconfig/lib" )
set( LIBCONFIG_LIBRARY "${prefix}/libconfig/src/libconfig-build/out/libconfig.so" )
set( LIBCONFIGPP_INCLUDE_DIR "${prefix}/libconfig/src/libconfig/lib" )
set( LIBCONFIGPP_LIBRARY "${prefix}/libconfig/src/libconfig-build/out/libconfig++.so" )
find_package(LibConfig REQUIRED)
list( APPEND dependencies libconfig )
endif()

find_package(NetCDF COMPONENTS CXX)
if( NOT NETCDF_FOUND )
message( "Ahoy, Captain! Downloading netcdf-c and netcdf-cxx4..." )
ExternalProject_Add( netcdf-c
                     PREFIX "${prefix}/netcdf-c"
                     GIT_REPOSITORY "git@github.com:Unidata/netcdf-c.git"
                     CMAKE_ARGS "-DENABLE_DAP=OFF" )

ExternalProject_Add( netcdf-cxx4
                     PREFIX "${prefix}/netcdf-cxx4"
                     GIT_REPOSITORY "git@github.com:Unidata/netcdf-cxx4.git" )

add_dependencies( netcdf-cxx4 netcdf-c )
set( NETCDF_LIBRARY "/usr/local/lib/libnetcdf.so" )
set( NETCDF_CXX_INCLUDE_DIR "${prefix}/netcdf-cxx4/src/netcdf-cxx4/cxx4" )
set( NETCDF_CXX_LIBRARY "/usr/local/lib/libnetcdf-cxx4.so" )
set( NETCDF_INCLUDE_DIR "${prefix}/netcdf-c/src/netcdf-c/include" )
find_package(NetCDF COMPONENTS CXX REQUIRED)
list( APPEND dependencies netcdf-c netcdf-cxx4 )
endif()
