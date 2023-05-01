# catch2

set( catch2_library_file_0
     "${prefix}/catch2-install/lib/libCatch2Main.a"
     CACHE FILEPATH ""
     FORCE )

set( catch2_library_file_1
     "${prefix}/catch2-install/lib/libCatch2.a"
     CACHE FILEPATH ""
     FORCE )

set( catch2_url "https://github.com/catchorg/Catch2.git" )

set( download_command git clone ${catch2_url} ${prefix}/catch2 )

set( configure_command
     ${CMAKE_COMMAND}
     -S ${prefix}/catch2
     -B ${prefix}/catch2-build
     -DCMAKE_INSTALL_PREFIX=${prefix}/catch2-install )

ExternalProject_Add( catch2_download
                     PREFIX ${prefix}
                     DOWNLOAD_COMMAND ${download_command}
                     CONFIGURE_COMMAND ${configure_command}
                     BUILD_BYPRODUCTS ${catch2_library_file_0} ${catch2_library_file_1}
                     BUILD_COMMAND ${CMAKE_COMMAND} 
                     --build ${prefix}/catch2-build --target install -- -j
                     INSTALL_COMMAND ${CMAKE_COMMAND} --install ${prefix}/catch2-build )

add_library( catch2 INTERFACE )

add_dependencies( catch2 catch2_download )

target_include_directories( catch2 INTERFACE ${prefix}/catch2-install/include )

target_link_libraries( catch2 INTERFACE 
                       ${catch2_library_file_0}
                       ${catch2_library_file_1} )

list( APPEND dependencies catch2 )
