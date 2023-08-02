set( cli_include_dir ${prefix}/cli11/include )

if( NOT EXISTS ${cli_include_dir} )

  set( cli11_url "https://github.com/CLIUtils/CLI11.git" )

  set( download_command 
       git clone --depth 1 --branch v2.1.2
       ${cli11_url} ${prefix}/cli11 )

  set( configure_command
      ${CMAKE_COMMAND} 
      -S ${prefix}/cli11
      -B ${prefix}/cli11_build
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_INSTALL_PREFIX=${prefix}/cli11_install )

  ExternalProject_Add( cli11_download
                       PREFIX ${prefix}
                       DOWNLOAD_COMMAND ${download_command}
                       CONFIGURE_COMMAND ""
                       BUILD_COMMAND ""
                       INSTALL_COMMAND "" )

endif()

# create an interface target so you can add a dependency to enforce_build_order
add_library( cli11 INTERFACE )

if( TARGET cli11_download )

  add_dependencies( cli11 cli11_download )

endif()

target_include_directories( cli11 INTERFACE ${cli_include_dir} )

list( APPEND dependencies cli11 )
