# Add thrust

if( NOT USE_CUDA )

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

add_library( thrust INTERFACE )

if( TARGET thrust_download )
  add_dependencies( thrust thrust_download )
endif()

include_directories( ${THRUST_INCLUDE_DIR} )

target_include_directories( thrust INTERFACE ${THRUST_INCLUDE_DIR} )
