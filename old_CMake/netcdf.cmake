# netcdf-c

set( NETCDF_INCLUDE_DIR "${prefix}/netcdf-c/include" CACHE PATH "" FORCE )

set( NETCDF_LIBRARY
     "${prefix}/netcdf-c-install/lib/libnetcdf${suffix}"
     CACHE FILEPATH ""
     FORCE )

set( NETCDF_CXX_INCLUDE_DIR "${prefix}/netcdf-cxx4-install/include" CACHE PATH "" FORCE)

set( NETCDF_CXX_LIBRARY
     "${prefix}/netcdf-cxx4-install/lib/libnetcdf-cxx4${suffix}"
     CACHE FILEPATH ""
     FORCE )

if ( NOT EXISTS ${NETCDF_INCLUDE_DIR} OR
     NOT EXISTS ${NETCDF_LIBRARY} OR
     NOT EXISTS ${NETCDF_CXX_INCLUDE_DIR} OR
     NOT EXISTS ${NETCDF_CXX_LIBRARY} )

  set( netcdf-c-url "https://github.com/Unidata/netcdf-c.git" )

  if( EXISTS ${prefix}/netcdf-c )

    set( download_command "" )

  else()

    set( download_command 
         git clone --depth 1 --branch v4.8.1
         ${netcdf-c-url}
         ${prefix}/netcdf-c )

    message( "netcdf-c will be downloaded..." )

  endif()

  set( configure_command
      ${CMAKE_COMMAND} 
      -S ${prefix}/netcdf-c
      -B ${prefix}/netcdf-c-build
      -DENABLE_DAP=OFF
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-c-install )

  ExternalProject_Add( netcdf-c_download
                       PREFIX ${prefix}
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ${configure_command}
                       BUILD_BYPRODUCTS ${NETCDF_LIBRARY}
                       BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-c-build -- -j
                       INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-c-build )

  # netcdf-cxx-4
  set( netcdf-cxx4-url "https://github.com/Unidata/netcdf-cxx4.git" )

  if( EXISTS ${prefix}/netcdf-cxx4 )

    set( download_command "" )

  else()

    set( download_command 
         git clone --depth 1 --branch v4.3.0
         ${netcdf-cxx4-url}
         ${prefix}/netcdf-cxx4 )

  endif()

  set( HDF5_PATHS ${HDF5_INCLUDE_DIRS} )

  string( REPLACE ";" ":" HDF5_PATHS "${HDF5_PATHS}" )

  set( configure_command
       PATH=${prefix}/netcdf-c-install/bin:${HDF5_PATHS}:$ENV{PATH}
       ${CMAKE_COMMAND}
       -S ${prefix}/netcdf-cxx4
       -B ${prefix}/netcdf-cxx4-build
       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
       -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-cxx4-install )

  ExternalProject_Add( netcdf-cxx4_download
                       PREFIX ${prefix}
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ${configure_command}
                       BUILD_BYPRODUCTS ${NETCDF_CXX_LIBRARY}
                       BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-cxx4-build -- -j
                       INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-cxx4-build )

  add_dependencies( netcdf-cxx4_download netcdf-c_download )

endif()

add_library( netcdf INTERFACE )

if( TARGET netcdf-c_download AND TARGET netcdf-cxx4_download )

  add_dependencies( netcdf netcdf-c_download )
  add_dependencies( netcdf netcdf-cxx4_download )

endif()

include_directories( ${NETCDF_CXX_INCLUDE_DIR} )

include_directories( ${NETCDF_INCLUDE_DIR} )

target_include_directories( netcdf INTERFACE
                            ${NETCDF_INCLUDE_DIR}
                            ${NETCDF_CXX_INCLUDE_DIR}
                            ${HDF5_INCLUDE_DIRS} )

target_link_libraries( netcdf INTERFACE
                       ${NETCDF_LIBRARY}
                       ${NETCDF_CXX_LIBRARY}
                       ${HDF5_LIBRARIES} )

list( APPEND dependencies netcdf )
