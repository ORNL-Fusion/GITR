
add_library( hdf5 INTERFACE )

# Captain! Find file:
find_file( hdf5_cxx_hl_shared_lib
           NAMES libhdf5_hl_cpp.so
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_file( hdf5_cxx_shared_lib
           NAMES libhdf5_cpp.so
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_file( hdf5_hl_shared_lib
           NAMES libhdf5_serial_hl.so
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_file( hdf5_shared_lib
           NAMES libhdf5_serial.so
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )



find_path( hdf5_headers
           NAMES hdf5.h
           HINTS "/usr/include/hdf5/serial" )

message( "Captain! hdf5_hl_shared_lib: ${hdf5_hl_shared_lib}" )

target_include_directories( hdf5 INTERFACE ${hdf5_headers} )

target_link_libraries( hdf5 INTERFACE 
                       ${hdf5_cxx_hl_shared_lib} 
                       ${hdf5_cxx_shared_lib}
                       ${hdf5_hl_shared_lib}
                       ${hdf5_shared_lib} )
