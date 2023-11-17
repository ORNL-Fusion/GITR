cmake -S $1 -B $2 -G Ninja \
-D NETCDF_C_HEADERS_DIR="/home/5n4/workshop/tree/external/netcdf-c-install/include" \
-D NETCDF_CXX_HEADERS_DIR="/home/5n4/workshop/tree/external/netcdf-cxx4-install/include" \
-D NETCDF_C_SHARED_LIB_DIR="/home/5n4/workshop/tree/external/netcdf-c-install/lib" \
-D NETCDF_CXX_SHARED_LIB_DIR="/home/5n4/workshop/tree/external/netcdf-cxx4-install/lib" \
-D LIBCONFIG_C_HEADERS_DIR="/home/5n4/workshop/tree/external/libconfig_install/include" \
-D LIBCONFIG_CXX_HEADERS_DIR="/home/5n4/workshop/tree/external/libconfig_install/include" \
-D LIBCONFIG_C_LIB_DIR="/home/5n4/workshop/tree/external/libconfig_install/lib" \
-D LIBCONFIG_CXX_LIB_DIR="/home/5n4/workshop/tree/external/libconfig_install/lib" \
&> $2/cmake_output.txt ;

cmake --build $2 -- -j 0 &> $2/build_output.txt
