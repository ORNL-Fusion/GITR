# Handle external dependencies

include( ExternalProject )

set( dependencies "" )

# Thrust

find_package(Thrust)

if( NOT THRUST_FOUND )

  if( NOT USE_CUDA )

    set( prefix "${CMAKE_CURRENT_SOURCE_DIR}/external" )

    message( "Downloading thrust..." )

    ExternalProject_Add( thrust_download
                         PREFIX "${prefix}/thrust" 
                         GIT_REPOSITORY "https://github.com/NVIDIA/thrust.git"
                         CONFIGURE_COMMAND ""
                         BUILD_COMMAND ""
                         INSTALL_COMMAND "" )

    set( THRUST_INCLUDE_DIR "${prefix}/thrust/src/thrust_download" CACHE PATH "" FORCE )

  else()

    set( THRUST_INCLUDE_DIR "${CUDAToolkit_INCLUDE_DIRS}" CACHE PATH "" FORCE )

  endif()

  set(THRUST_INCLUDE_DIRS ${THRUST_INCLUDE_DIR} CACHE PATH "" FORCE )

  find_package(Thrust REQUIRED)

endif()

add_library( thrust INTERFACE )

if( TARGET thrust_download )
  add_dependencies( thrust thrust_download )
endif()

include_directories( ${THRUST_INCLUDE_DIR} )

target_include_directories( thrust INTERFACE ${THRUST_INCLUDE_DIR} )

list( APPEND dependencies thrust )

# Libconfig

find_package(LibConfig)

if( NOT LIBCONFIG_FOUND )

  message( "Downloading libconfig..." )

  ExternalProject_Add( libconfig_download
                       PREFIX "${prefix}/libconfig"
                       GIT_REPOSITORY "https://github.com/hyperrealm/libconfig.git" 
                       INSTALL_COMMAND "" )

  set( LIBCONFIG_INCLUDE_DIR 
       "${prefix}/libconfig/src/libconfig_download/lib" 
       CACHE PATH "" FORCE )


  set( LIBCONFIG_LIBRARY 
       "${prefix}/libconfig/src/libconfig_download-build/out/libconfig.so" 
       CACHE FILEPATH "" FORCE )

  set( LIBCONFIGPP_INCLUDE_DIR 
       "${prefix}/libconfig/src/libconfig_download/lib"
       CACHE PATH "" FORCE )


  set( LIBCONFIGPP_LIBRARY
       "${prefix}/libconfig/src/libconfig_download-build/out/libconfig++.so"
       CACHE FILEPATH "" FORCE )

  find_package(LibConfig REQUIRED)

endif()

add_library( libconfig INTERFACE )

if( TARGET libconfig_download )
  add_dependencies( libconfig libconfig_download )
endif()

include_directories( ${LIBCONFIG_INCLUDE_DIR} )

include_directories( ${LIBCONFIGPP_INCLUDE_DIR} )

target_include_directories( libconfig INTERFACE 
                            ${LIBCONFIG_INCLUDE_DIR}
                            ${LIBCONFIGPP_INCLUDE_DIR} )

target_link_libraries( libconfig INTERFACE
                       ${LIBCONFIG_LIBRARY}
                       ${LIBCONFIGPP_LIBRARY} )

list( APPEND dependencies libconfig )

# NETCDF

find_package(NetCDF COMPONENTS CXX)

if( NOT NETCDF_FOUND )

  message( "Downloading netcdf-c and netcdf-cxx4..." )

  ExternalProject_Add( netcdf-c_download
                       PREFIX "${prefix}/netcdf-c"
                       GIT_REPOSITORY "https://github.com/Unidata/netcdf-c.git"
                       CMAKE_ARGS "-DENABLE_DAP=OFF" )

  ExternalProject_Add( netcdf-cxx4_download
                       PREFIX "${prefix}/netcdf-cxx4"
                       GIT_REPOSITORY "https://github.com/Unidata/netcdf-cxx4.git" )

  add_dependencies( netcdf-cxx4_download netcdf-c_download )

  set( NETCDF_LIBRARY "/usr/local/lib/libnetcdf.so" CACHE FILEPATH "" FORCE )

  set( NETCDF_CXX_INCLUDE_DIR "${prefix}/netcdf-cxx4/src/netcdf-cxx4_download/cxx4" CACHE PATH "" FORCE)

  set( NETCDF_CXX_LIBRARY "/usr/local/lib/libnetcdf-cxx4.so" CACHE FILEPATH "" FORCE )

  set( NETCDF_INCLUDE_DIR "${prefix}/netcdf-c/src/netcdf-c_download/include" CACHE PATH "" FORCE )

  find_package(NetCDF COMPONENTS CXX REQUIRED)

endif()

add_library( netcdf INTERFACE )

if( TARGET netcdf-c_download AND TARGET netcdf-cxx4_download )
  add_dependencies( netcdf netcdf_c_download netcdf-cxx4_download )
endif()

include_directories( ${NETCDF_CXX_INCLUDE_DIR} )

include_directories( ${NETCDF_INCLUDE_DIR} )

target_include_directories( netcdf INTERFACE ${NETCDF_INCLUDE_DIR} ${NETCDF_CXX_INCLUDE_DIR} )

target_link_libraries( netcdf INTERFACE ${NETCDF_LIBRARY} ${NETCDF_CXX_LIBRARY} )

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
