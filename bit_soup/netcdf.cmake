# Add an interface library here for netcdf

add_library( netcdf INTERFACE )

target_include_directories( netcdf INTERFACE "/usr/include" )

target_link_libraries( netcdf INTERFACE "/usr/lib/libnetcdf_c++4.so" "/usr/lib/libnetcdf.so" )
