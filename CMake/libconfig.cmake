add_library( libconfig INTERFACE )

find_file( libconfig_cxx_lib
           NAMES libconfig++.so 
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_file( libconfig_c_lib 
           NAMES libconfig.so
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_path( libconfig_headers 
           NAMES libconfig.h )

message( "Ahoy! libconfig_cxx_lib: ${libconfig_cxx_lib}" )
message( "Ahoy! libconfig_c_lib: ${libconfig_c_lib}" )

target_include_directories( libconfig INTERFACE ${libconfig_headers} )

target_link_libraries( libconfig INTERFACE "${libconfig_cxx_lib}" "${libconfig_c_lib}" )
