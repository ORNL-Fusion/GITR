# Add thrust

if( NOT GITR_USE_CUDA )

  set( THRUST_INCLUDE_DIR "${prefix}/thrust" CACHE PATH "" FORCE )

  if( NOT EXISTS ${THRUST_INCLUDE_DIR} )
    
    message( "thrust will be downloaded..." )

    set( thrust_url "https://github.com/NVIDIA/thrust.git" )

    set( download_command
         git clone --depth 1 --branch 1.15.0
         ${thrust_url} 
         ${prefix}/thrust )

    ExternalProject_Add( thrust_download
                         PREFIX ${prefix} 
                         DOWNLOAD_COMMAND ${download_command}
                         CONFIGURE_COMMAND ""
                         BUILD_COMMAND ""
                         INSTALL_COMMAND "" )

  endif()

else()

  set( THRUST_INCLUDE_DIR "${CUDAToolkit_INCLUDE_DIRS}" CACHE PATH "" FORCE )

endif()

set(THRUST_INCLUDE_DIRS ${THRUST_INCLUDE_DIR} CACHE PATH "" FORCE )

add_library( thrust INTERFACE )

if( TARGET thrust_download )

  add_dependencies( thrust thrust_download )

endif()

include_directories( ${THRUST_INCLUDE_DIR} )

target_include_directories( thrust INTERFACE ${THRUST_INCLUDE_DIR} )
