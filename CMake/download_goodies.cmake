include( ExternalProject )
# functions to obtain 3rd party dependencies
function( download_goodies )

# no prefix is needed in this simple function
# obtain arguments to pass to the ExternalProject_Add
# an empty install command obviates the install step
set( prepend_string "external" )
set( bool_options "" )
#Captain! Can you move this into a list append? Or do it with multiple strings? I bet
set( key_single_val "name;prefix;git_repository;install_command;build_command" )
set( key_multi_val "cmake_args" )

cmake_parse_arguments( "${prepend_string}"
                       "${bool_options}"
                       "${key_single_val}"
                       "${key_multi_val}"
                       "${ARGN}" )

message( STATUS "build_command: ${external_build_command}")
message( STATUS "install_command: ${external_install_command}")

if( NOT DEFINED external_build_command )
  message( STATUS "Captain! external_build_command not defined. Setting to empty string:")
  set( external_build_command "" )
  if( DEFINED external_build_command )
    message( STATUS "external_build_command set to: ${external_build_command}")
  endif()
endif()

if( NOT DEFINED external_install_command )
  message( STATUS "Captain! external_install_command not defined. Setting to empty string:")
  set( external_install_command "" )
  if( DEFINED external_install_command )
    message( STATUS "external_install_command set to: ${external_install_command}")
  endif()
endif()

# empty strings handled differently than variables set to empty strings
ExternalProject_Add( ${external_name} 
                     PREFIX ${external_prefix}
                     GIT_REPOSITORY ${external_git_repository}
                     CMAKE_ARGS ${external_cmake_args}
                     BUILD_COMMAND ${external_build_command}
                     INSTALL_COMMAND ${external_install_command} )

endfunction()
