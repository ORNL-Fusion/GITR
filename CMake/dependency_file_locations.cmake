
find_file( netcdf_cxx_shared_lib 
           NAMES libnetcdf_c++4.so 
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_file( netcdf_c_shared_lib 
           NAMES libnetcdf.so 
           HINTS "/usr/lib/x86_64-linux-gnu" "/usr/lib" )

find_path( netcdf_headers
           NAMES netcdf
           HINTS "/usr/include" )

# stuff for libconfig etc

# user just modifies this file - all dependencies must be filled in

