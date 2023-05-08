cmake -S $1 -B $2 -G Ninja \
-D NETCDF_HEADERS_DIR="/host/GITR" \
-D NETCDF_CXX_SHARED_LIB_DIR="/host/tree/build" \ 
-D NETCDF_C_SHARED_LIB_DIR="/host/tree/build" \ 
-D LIBCONFIG_CXX_LIB_DIR="/host/tree/build" \ 
-D NETCDF_C_LIB_DIR="/host/tree/build" \ 
&> $2/cmake_output.txt ;

cmake --build $2 -- -j 2 &> $2/build_output.txt
