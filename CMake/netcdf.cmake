# netcdf-c

find_package(NetCDF COMPONENTS CXX)

if( NOT NETCDF_FOUND )

  message( "Downloading netcdf-c and netcdf-cxx4..." )

  set( netcdf-c-url "https://github.com/Unidata/netcdf-c.git" )

  # Captain! Set the HDF5_ROOT or other hdf5 variables so netcdf will link to it
  set( configure_command
      ${CMAKE_COMMAND} 
      -S ${prefix}/netcdf-c
      -B ${prefix}/netcdf-c-build
      -DENABLE_DAP=OFF
      -DHDF5_ROOT=${hdf5_install_dir}
      -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-c-install
      -DCMAKE_INSTALL_RPATH=${prefix}/netcdf-c-install/lib )

  if( EXISTS ${prefix}/netcdf-c )

    set( download_command "" )

  else()

    set( download_command git clone ${netcdf-c-url} ${prefix}/netcdf-c )

  endif()

  ExternalProject_Add( netcdf-c_download
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ${configure_command}
                       BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-c-build -- -j
                       INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-c-build )

  # netcdf-cxx-4

  set( netcdf-cxx4-url "https://github.com/Unidata/netcdf-cxx4.git" )

  set( configure_command
       PATH=${prefix}/netcdf-c-install/bin:$ENV{PATH}
       ${CMAKE_COMMAND}
       -S ${prefix}/netcdf-cxx4
       -B ${prefix}/netcdf-cxx4-build
       -DCMAKE_INSTALL_PREFIX=${prefix}/netcdf-cxx4-install )

  if( EXISTS ${prefix}/netcdf-cxx4 )

    set( download_command "" )

  else()

    set( download_command git clone ${netcdf-cxx4-url} ${prefix}/netcdf-cxx4 )

  endif()

  ExternalProject_Add( netcdf-cxx4_download
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ${configure_command}
                       BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/netcdf-cxx4-build -- -j
                       INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/netcdf-cxx4-build )

  add_dependencies( netcdf-cxx4_download netcdf-c_download )

  set( NETCDF_LIBRARY
       "${prefix}/netcdf-c-install/lib/libnetcdf${suffix}"
       CACHE FILEPATH ""
       FORCE )

  set( NETCDF_CXX_INCLUDE_DIR "${prefix}/netcdf-cxx4-install/include" CACHE PATH "" FORCE)

  set( NETCDF_CXX_LIBRARY
       "${prefix}/netcdf-cxx4-install/lib/libnetcdf-cxx4${suffix}"
       CACHE FILEPATH ""
       FORCE )

  set( NETCDF_INCLUDE_DIR "${prefix}/netcdf-c/include" CACHE PATH "" FORCE )

  find_package(NetCDF COMPONENTS CXX REQUIRED)

endif()

add_library( netcdf INTERFACE )

if( TARGET netcdf-c_download AND TARGET netcdf-cxx4_download )

  add_dependencies( netcdf-c_download hdf5 )
  add_dependencies( netcdf-cxx4_download hdf5 )
  add_dependencies( netcdf netcdf_c_download netcdf-cxx4_download hdf5 )

endif()

# Captain! Are these global includes really necessary? Comment out to test
# include_directories( ${NETCDF_CXX_INCLUDE_DIR} )

# include_directories( ${NETCDF_INCLUDE_DIR} )

target_include_directories( netcdf INTERFACE
                            ${NETCDF_INCLUDE_DIR}
                            ${NETCDF_CXX_INCLUDE_DIR} )
                            # Captain! Will the target link libraries call link correctly?
                            #${HDF5_INCLUDE_DIRS} )

target_link_libraries( netcdf INTERFACE
                       ${NETCDF_LIBRARY}
                       ${NETCDF_CXX_LIBRARY}
                       hdf5 )
                       # Will the target handle this correctly?
                       #${HDF5_LIBRARIES} )
