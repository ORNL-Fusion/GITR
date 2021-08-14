# Add thrust

if( NOT GITR_USE_CUDA )

  message( "Downloading thrust..." )

  set( thrust_url "https://github.com/ORNL-Fusion/thrust_archive.git" )

  set( download_command git clone ${thrust_url} ${prefix}/thrust )

  ExternalProject_Add( thrust_download
                       PREFIX ${prefix} 
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ""
                       BUILD_COMMAND ""
                       INSTALL_COMMAND "" )

  set( THRUST_INCLUDE_DIR "${prefix}/thrust" CACHE PATH "" FORCE )

else()

  set( THRUST_INCLUDE_DIR "${CUDAToolkit_INCLUDE_DIRS}" CACHE PATH "" FORCE )

endif()

set(THRUST_INCLUDE_DIRS ${THRUST_INCLUDE_DIR} CACHE PATH "" FORCE )

find_package(Thrust REQUIRED)

add_library( thrust INTERFACE )

if( TARGET thrust_download )

  add_dependencies( thrust thrust_download )

endif()

include_directories( ${THRUST_INCLUDE_DIR} )

target_include_directories( thrust INTERFACE ${THRUST_INCLUDE_DIR} )
