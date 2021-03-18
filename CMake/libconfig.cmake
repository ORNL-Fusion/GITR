# Add libconfig

find_package(LibConfig)

if( NOT LIBCONFIG_FOUND )

  message( "Downloading libconfig..." )

  set( libconfig_url "https://github.com/hyperrealm/libconfig.git" )

  if( EXISTS ${prefix}/libconfig )
    ExternalProject_Add( libconfig_download
                         DOWNLOAD_COMMAND ""
                         CONFIGURE_COMMAND cmake -S ${prefix}/libconfig -B ${prefix}/libconfig_build -DCMAKE_INSTALL_PREFIX=${prefix}/libconfig_install
                         BUILD_COMMAND cmake --build ${prefix}/libconfig_build -- -j
                         INSTALL_COMMAND cmake --install ${prefix}/libconfig_build ) 
  else()
    ExternalProject_Add( libconfig_download
                         DOWNLOAD_COMMAND git clone ${libconfig_url} ${prefix}/libconfig
                         CONFIGURE_COMMAND cmake -S ${prefix}/libconfig -B ${prefix}/libconfig_build -DCMAKE_INSTALL_PREFIX=${prefix}/libconfig_install
                         BUILD_COMMAND cmake --build ${prefix}/libconfig_build -- -j
                         INSTALL_COMMAND cmake --install ${prefix}/libconfig_build ) 
  endif()


  set( LIBCONFIG_INCLUDE_DIR 
       "${prefix}/libconfig_install/include" 
       CACHE PATH "" FORCE )


  set( LIBCONFIG_LIBRARY 
       "${prefix}/libconfig_install/lib/libconfig.so" 
       CACHE FILEPATH "" FORCE )

  set( LIBCONFIGPP_INCLUDE_DIR 
       "${prefix}/libconfig_install/include" 
       CACHE PATH "" FORCE )


  set( LIBCONFIGPP_LIBRARY
       "${prefix}/libconfig_install/lib/libconfig++.so" 
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
