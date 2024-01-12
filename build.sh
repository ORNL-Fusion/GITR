cmake -S $1 -B $2 -G Ninja \
-DNETCDF_C_HEADERS_DIR="/home/5n4/build_right/external/netcdf-c-install/include" \
-DNETCDF_CXX_HEADERS_DIR="/home/5n4/build_right/external/netcdf-cxx4-install/include" \
-DNETCDF_C_SHARED_LIB_DIR="/home/5n4/build_right/external/netcdf-c-install/lib" \
-DNETCDF_CXX_SHARED_LIB_DIR="/home/5n4/build_right/external/netcdf-cxx4-install/lib" \
-DLIBCONFIG_C_HEADERS_DIR="/home/5n4/build_right/external/libconfig_install/include" \
-DLIBCONFIG_CXX_HEADERS_DIR="/home/5n4/build_right/external/libconfig_install/include" \
-DLIBCONFIG_C_LIB_DIR="/home/5n4/build_right/external/libconfig_install/lib" \
-DLIBCONFIG_CXX_LIB_DIR="/home/5n4/build_right/external/libconfig_install/lib" \
&> $2/cmake_output.txt ;

cmake --build $2 -- -j 0 &> $2/build_output.txt
