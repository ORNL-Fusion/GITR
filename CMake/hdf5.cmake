
# check if file already exists:
if( NOT EXISTS ${prefix}/CMake-hdf5-1.13.0/build/bin/libhdf5${suffix} AND 
    NOT EXISTS ${prefix}/CMake-hdf5-1.13.0/hdf5-1.13.0/src/hdf5.h ) 

set( hdf5_archive_url
     "https://raw.githubusercontent.com/ORNL-Fusion/hdf5_archive/main/CMake-hdf5-1.13.0.tar.gz" )

set( hdf5_archive ${prefix}/hdf5.tar.gz )

if( NOT EXISTS ${hdf5_archive} )


  file( DOWNLOAD ${hdf5_archive_url} ${hdf5_archive} )

endif()

file( ARCHIVE_EXTRACT 
      INPUT ${hdf5_archive}
      DESTINATION ${prefix} )

add_custom_target( hdf5_build )

set( pre_build_command 
     ${CMAKE_CTEST_COMMAND}
     -S HDF5config.cmake,BUILD_GENERATOR=Unix
     -C Release
     -V
     -O
     hdf5.log )

add_custom_command( TARGET hdf5_build
                    PRE_BUILD COMMAND ${pre_build_command}
                    WORKING_DIRECTORY ${prefix}/CMake-hdf5-1.13.0 )

endif()

add_library( hdf5 INTERFACE )

# Captain! Check for hdf5_build here before adding it as a dependency
add_dependencies( hdf5 hdf5_build )

target_include_directories( hdf5 
                            INTERFACE 
                            ${prefix}/CMake-hdf5-1.13.0/hdf5-1.13.0/src 
                            ${prefix}/CMake-hdf5-1.13.0/build/src )

target_link_libraries( hdf5 INTERFACE ${prefix}/CMake-hdf5-1.13.0/build/bin/libhdf5${suffix} )
