# Add an interface library here for netcdf

add_library( netcdf INTERFACE )

# Captain! Find file:
find_file( netcdf_cxx_shared_lib 
           NAMES libnetcdf_c++4.so 
           HINTS ${NETCDF_CXX_SHARED_LIB_DIR} )

find_file( netcdf_c_shared_lib 
           NAMES libnetcdf.so 
           HINTS ${NETCDF_C_SHARED_LIB_DIR} )

find_path( netcdf_headers
           NAMES netcdf
           HINTS ${NETCDF_HEADERS_DIR} )

target_include_directories( netcdf INTERFACE ${netcdf_headers} )

target_link_libraries( netcdf INTERFACE ${netcdf_cxx_shared_lib} ${netcdf_c_shared_lib} )
