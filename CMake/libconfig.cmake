add_library( libconfig INTERFACE )

# Captain! It is untested whether this works in alpine too
find_file( libconfig_cxx_lib
           NAMES libconfig++.so 
           HINTS ${LIBCONFIG_CXX_LIB_DIR} )

find_file( libconfig_c_lib 
           NAMES libconfig.so
           HINTS ${LIBCONFIG_C_LIB_DIR} )

find_path( libconfig_c_headers 
           NAMES libconfig.h 
           HINTS ${LIBCONFIG_C_HEADERS_DIR} )

find_path( libconfig_cxx_headers 
           NAMES libconfig.h++
           HINTS ${LIBCONFIG_CXX_HEADERS_DIR} )

target_include_directories( libconfig INTERFACE ${libconfig_c_headers} ${libconfig_cxx_headers} )

target_link_libraries( libconfig INTERFACE ${libconfig_cxx_lib} ${libconfig_c_lib} )
