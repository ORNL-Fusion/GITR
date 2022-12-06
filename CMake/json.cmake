# Captain! Add Ninja as the default cmake generator for any of these projects
set( JSON_BuildTests OFF CACHE INTERNAL "" )

set( json_url "https://github.com/nlohmann/json.git" )

set( json_file ${prefix}/json/json_install/include/nlohmann/json.hpp )

set( download_command 
     git clone --depth 1 --branch v3.10.5
     ${json_url}
     ${prefix}/json )

set( configure_command
    ${CMAKE_COMMAND} 
    -S ${prefix}/json
    -B ${prefix}/json_build
    -DENABLE_DAP=OFF
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_INSTALL_PREFIX=${prefix}/json_install )

# Download only
ExternalProject_Add( json_download
                     PREFIX ${prefix}
                     DOWNLOAD_COMMAND ${download_command}
                     CONFIGURE_COMMAND ${configure_command}
                     BUILD_BYPRODUCTS ${json_file}
                     BUILD_COMMAND ${CMAKE_COMMAND} --build ${prefix}/json_build -- -j
                     INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/json_build )



add_library( json_lib INTERFACE )

if( TARGET json_download )

  add_dependencies( json_lib json_download )

endif()

target_include_directories( json_lib INTERFACE ${prefix}/json_install/include )

# Captain!
add_executable( json src/json.cpp )

add_dependencies( json json_lib )

target_link_libraries( json json_lib )
