# find files

#message( STATUS "libconfig_cxx_lib: libconfig++.${suffix}" )

find_file( libconfig_cxx_lib
           NAMES libconfig++${suffix}
           HINTS ${LIBCONFIG_CXX_LIB_DIR} )

find_file( libconfig_c_lib 
           NAMES libconfig${suffix}
           HINTS ${LIBCONFIG_C_LIB_DIR} )

find_path( libconfig_c_headers 
           NAMES libconfig.h 
           HINTS ${LIBCONFIG_C_HEADERS_DIR} )

find_path( libconfig_cxx_headers 
           NAMES libconfig.h++
           HINTS ${LIBCONFIG_CXX_HEADERS_DIR} )

# create imported target out of the files
add_library( libconfig_cxx SHARED IMPORTED )
add_library( libconfig_c SHARED IMPORTED )

set_property( TARGET libconfig_cxx PROPERTY IMPORTED_LOCATION ${libconfig_cxx_lib} )
set_property( TARGET libconfig_c PROPERTY IMPORTED_LOCATION ${libconfig_c_lib} )

target_include_directories( libconfig_cxx 
                            INTERFACE 
                            ${libconfig_cxx_headers} )

target_include_directories( libconfig_c 
                            INTERFACE 
                            ${libconfig_c_headers} )

include_directories( ${libconfig_cxx_headers} ${libconfig_c_headers} )
