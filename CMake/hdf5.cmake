# Captain! Replace the externalproject with a file download etc - externalproject is overkill
# in this situation
# Download and extract hdf5 archived file
set( hdf5_archive_url
    "https://github.com/ORNL-Fusion/hdf5_archive/blob/main/CMake-hdf5-1.12.1.tar.gz?raw=true" )

ExternalProject_Add( hdf5_download 
                     PREFIX ${prefix}/hdf5
                     URL ${hdf5_archive_url}
                     CONFIGURE_COMMAND ""
                     BUILD_COMMAND  ""
                     INSTALL_COMMAND "" )

# Build the tarball full of the build output
cmake_path( GET CMAKE_CTEST_COMMAND PARENT_PATH parent_path )

# Build szlibH
# set( szlib_directory ${prefix}/hdf5/src/)
# add_custom_command( OUTPUT
#                     COMMAND
#                     WORKING_DIRECTORY )

set( tarball ${prefix}/hdf5/src/hdf5_download/HDF5-1.12.1-Linux.tar.gz )

set( command PATH=$ENV{PATH}:${parent_path} bash ${prefix}/hdf5/src/hdf5_download/build-unix.sh )

set( working_directory ${prefix}/hdf5/src/hdf5_download )

add_custom_command( OUTPUT ${tarball}
                    COMMAND ${command}
                    WORKING_DIRECTORY ${working_directory} )

add_custom_target( hdf5_build DEPENDS ${tarball} )

add_dependencies( hdf5_build hdf5_download )

# Captain! Can you now just turn this into a custom target as before? Without referencing files?
set( build_artifacts ${prefix}/hdf5/src/hdf5_download/HDF5-1.12.1-Linux )

add_custom_command( OUTPUT ${build_artifacts}
                    COMMAND 
                    ${CMAKE_COMMAND} -E 
                    tar xzf ${prefix}/hdf5/src/hdf5_download/HDF5-1.12.1-Linux.tar.gz
                    WORKING_DIRECTORY ${prefix}/hdf5 )

add_custom_target( hdf5_extract DEPENDS ${build_artifacts} )

add_dependencies( hdf5_extract hdf5_build )

# Captain! None of this is used
set( hdf5_include_dirs
     ${prefix}/hdf5/HDF5-1.12.1-Linux/HDF_Group/HDF5/1.12.1/include )

set( hdf5_libraries
     ${prefix}/hdf5/HDF5-1.12.1-Linux/HDF_Group/HDF5/1.12.1/lib/libhdf5.so )

# Captain! You might just want to go ahead and hardcode the variables in FindHDF5 and pass
# them as -D options instead of HDF5_root in the netcdf build command
add_library( hdf5 INTERFACE )

# How do we ensure that it is downloaded before cmake attempts to build?
add_dependencies( hdf5 hdf5_build )

add_dependencies( hdf5 hdf5_extract )
