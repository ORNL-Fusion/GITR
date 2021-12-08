set( cli11_url "https://github.com/CLIUtils/CLI11.git" )

set( download_command git clone ${cli11_url} ${prefix}/cli11 )

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
                     CONFIGURE_COMMAND ${configure_command}
                     BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/cli11_build -- -j
                     INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/cli11_build )

# create an interface target so you can add a dependency to enforce_build_order
add_library( cli11 INTERFACE )

add_dependencies( cli11 cli11_download )

target_include_directories( cli11 INTERFACE ${prefix}/cli11_install/include )
