# Captain! Revert back to downloading the r
# Captain! This simply forces all the steps to run for testing/development purposes.
add_library( hdf5 INTERFACE )

# Obtain the archived repo with tarballs in it
set( hdf5_git_repo "https://github.com/ORNL-Fusion/hdf5_archive.git" )

ExternalProject_Add( hdf5_download
                     PREFIX ${prefix}/hdf5
                     GIT_REPOSITORY ${hdf5_git_repo}
                     GIT_TAG main
                     CONFIGURE_COMMAND ""
                     BUILD_COMMAND  ""
                     INSTALL_COMMAND "" )

# Build the szlib tarball
set( szlib_repo_tar ${prefix}/hdf5/src/hdf5_download/szip-2.1.1.tar.gz )
set( hdf5_repo_tar ${prefix}/hdf5/src/hdf5_download/CMake-hdf5-1.12.1.tar.gz )
set( zlib_repo_tar ${prefix}/hdf5/src/hdf5_download/zlib.tar.gz )
set( szlib_repo ${prefix}/hdf5/src/hdf5_download/szip-2.1.1 )
set( hdf5_repo ${prefix}/hdf5/src/hdf5_download/CMake-hdf5-1.12.1 )
set( zlib_repo ${prefix}/hdf5/src/hdf5_download/ZLib )

# Create a custom command to extract them
add_custom_command( OUTPUT ${szlib_repo} ${hdf5_repo} ${zlib_repo}
                    COMMAND ${CMAKE_COMMAND} -E 
                    tar -xzf ${szlib_repo_tar}
                    COMMAND ${CMAKE_COMMAND} -E 
                    tar -xzf ${hdf5_repo_tar}
                    COMMAND ${CMAKE_COMMAND} -E 
                    tar -xzf ${zlib_repo_tar}
                    WORKING_DIRECTORY ${prefix}/hdf5/src/hdf5_download )

add_custom_target( extract_repo_contents DEPENDS ${szlib_repo} ${hdf5_repo} ${zlib_repo} )

add_dependencies( extract_repo_contents hdf5_download )

add_dependencies( hdf5 extract_repo_contents )

# Create a custom target to build szlib
set( szip_install ${prefix}/szip )

add_custom_command( OUTPUT ${szip_install}
                    COMMAND ./configure --prefix=${szip_install}
                    COMMAND make
                    COMMAND make install
                    WORKING_DIRECTORY ${prefix}/hdf5/src/hdf5_download/szip-2.1.1 )

add_custom_target( build_szip DEPENDS ${szip_install} )

add_dependencies( build_szip extract_repo_contents )

add_dependencies( hdf5 build_szip )

# Create a custom target to build zlib...
set( zlib_install ${prefix}/zlib )

add_custom_command( OUTPUT ${zlib_install} 
                    COMMAND ./configure --prefix=${zlib_install}
                    COMMAND make
                    COMMAND make install
                    WORKING_DIRECTORY ${prefix}/hdf5/src/hdf5_download/ZLib )

add_custom_target( build_zlib DEPENDS ${zlib_install} )

add_dependencies( build_zlib extract_repo_contents )

add_dependencies( hdf5 build_zlib )

# Create a custom target to build hdf5 - fingers crossed this actually gets built.
set( hdf5_install_tar ${prefix}/hdf5/install/HDF5-1.12.1-Linux.tar.gz )

set( working_directory ${prefix}/hdf5/src/hdf5_download/CMake-hdf5-1.12.1 )

set( command ${CMAKE_CTEST_COMMAND}
     -DADD_BUILD_OPTIONS="-DSZIP_LIBRARY=${szip_install}/libsz.so -DSZIP_INCLUDE_DIR=${szip_install}/include -DZLIB_LIBRARY=${zlib_install}/lib/libz.so -DZLIB_INCLUDE_DIR=${zlib_install}/include" 
     -S HDF5config.cmake,INSTALLDIR=${prefix}/hdf5/install,BUILD_GENERATOR=Unix
     -C Release -V -O hdf5.log && cd ${working_directory}/build && make install )

# Captain! Test out whether you can add multiple "working directory" lines at the bottom - 
# simply adding another command to "make install" since the ctest one doesn't do it...
add_custom_command( OUTPUT ${hdf5_install_tar} 
                    COMMAND ${command} 
                    WORKING_DIRECTORY ${working_directory} )

set( hdf5_lib ${prefix}/hdf5/install/lib )
set( hdf5_include ${prefix}/hdf5/install/include )

add_custom_target( build_hdf5 DEPENDS ${hdf5_install_tar} )

add_dependencies( build_hdf5 build_szip )

add_dependencies( build_hdf5 build_zlib )

add_dependencies( hdf5 build_hdf5 )

# An ldd on the libhdf5.so should indicate that it is linked to the szip built here.
# You can then use the netcdf-c CMakeLists.txt:747 work by passing the right variable to CMake.
# with SZIP_LIBRARY. Link it to the same one and then see if you can get it to actually compile.





















