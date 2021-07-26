# Add the netcdf library

find_package(NetCDF COMPONENTS CXX)

if( NOT NETCDF_FOUND )

  message( "Downloading netcdf-c and netcdf-cxx4..." )

  set( netcdf-c-url "https://github.com/Unidata/netcdf-c.git" )

  if( EXISTS ${prefix}/netcdf-c )
    ExternalProject_Add( netcdf-c_download
                         DOWNLOAD_COMMAND ""
                         CONFIGURE_COMMAND ${CMAKE_COMMAND} -S ${prefix}/netcdf-c -B ${prefix}/netcdf-c-build -DENABLE_DAP=OFF -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-c-install
                         BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-c-build -- -j
                         INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-c-build )
  else()
    ExternalProject_Add( netcdf-c_download
                         DOWNLOAD_COMMAND git clone ${netcdf-c-url} ${prefix}/netcdf-c
                         CONFIGURE_COMMAND ${CMAKE_COMMAND} -S ${prefix}/netcdf-c -B ${prefix}/netcdf-c-build -DENABLE_DAP=OFF -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-c-install
                         BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-c-build -- -j
                         INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-c-build )

  endif()

  set( netcdf-cxx4-url "https://github.com/Unidata/netcdf-cxx4.git" )
  if( EXISTS ${prefix}/netcdf-cxx4 )
    # No need to install this one
    ExternalProject_Add( netcdf-cxx4_download
                         DOWNLOAD_COMMAND ""
                         CONFIGURE_COMMAND PATH=${prefix}/netcdf-c-install/bin:$ENV{PATH} ${CMAKE_COMMAND} -S ${prefix}/netcdf-cxx4 -B ${prefix}/netcdf-cxx4-build -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-cxx4-install
                         BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-cxx4-build -- -j
                         INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-cxx4-build )
  else()
    # No need to install this one
    ExternalProject_Add( netcdf-cxx4_download
                         DOWNLOAD_COMMAND git clone ${netcdf-cxx4-url} ${prefix}/netcdf-cxx4
                         CONFIGURE_COMMAND PATH=${prefix}/netcdf-c-install/bin:$ENV{PATH} ${CMAKE_COMMAND} -S ${prefix}/netcdf-cxx4 -B ${prefix}/netcdf-cxx4-build -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-cxx4-install
                         BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-cxx4-build -- -j
                         INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-cxx4-build )
  endif()

  add_dependencies( netcdf-cxx4_download netcdf-c_download )

  set( NETCDF_LIBRARY "${prefix}/netcdf-c-install/lib/libnetcdf${suffix}" CACHE FILEPATH "" FORCE )

  set( NETCDF_CXX_INCLUDE_DIR "${prefix}/netcdf-cxx4-install/include" CACHE PATH "" FORCE)

  set( NETCDF_CXX_LIBRARY "${prefix}/netcdf-cxx4-install/lib/libnetcdf-cxx4${suffix}" CACHE FILEPATH "" FORCE )

  set( NETCDF_INCLUDE_DIR "${prefix}/netcdf-c/include" CACHE PATH "" FORCE )

  find_package(NetCDF COMPONENTS CXX REQUIRED)

endif()

add_library( netcdf INTERFACE )

if( TARGET netcdf-c_download AND TARGET netcdf-cxx4_download )
  add_dependencies( netcdf netcdf_c_download netcdf-cxx4_download )
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


























