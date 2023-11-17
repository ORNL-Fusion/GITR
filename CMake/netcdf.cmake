# find files

# Captain! Find file:
find_file( netcdf_cxx_shared_lib 
              NAMES libnetcdf-cxx4.so libnetcdf_c++4.so
              HINTS ${NETCDF_CXX_SHARED_LIB_DIR} )

find_file( netcdf_c_shared_lib 
              NAMES libnetcdf.so 
              HINTS ${NETCDF_C_SHARED_LIB_DIR} )

find_path( netcdf_c_headers
           NAMES netcdf.h
           HINTS ${NETCDF_C_HEADERS_DIR} )

find_path( netcdf_cxx_headers
           NAMES netcdf
           HINTS ${NETCDF_CXX_HEADERS_DIR} )

message( "netcdf_c_shared_lib: ${netcdf_c_shared_lib}" )
add_library( netcdf_cxx SHARED IMPORTED )
add_library( netcdf_c SHARED IMPORTED )

# create imported library with the files
set_property( TARGET netcdf_cxx PROPERTY IMPORTED_LOCATION ${netcdf_cxx_shared_lib} )
set_property( TARGET netcdf_c PROPERTY IMPORTED_LOCATION ${netcdf_c_shared_lib} )
target_include_directories( netcdf_cxx INTERFACE ${netcdf_cxx_headers} )
target_include_directories( netcdf_c INTERFACE ${netcdf_c_headers} )


include_directories( ${netcdf_cxx_headers} ${netcdf_c_headers} )
