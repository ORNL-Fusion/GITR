cmake -S $1 -B $2 \
-DNETCDF_C_HEADERS_DIR="/opt/cray/pe/netcdf/4.9.0.9/nvidia/23.3/include" \
-DNETCDF_CXX_HEADERS_DIR="/opt/cray/pe/netcdf/4.9.0.9/nvidia/23.3/include" \
-DNETCDF_C_SHARED_LIB_DIR="/opt/cray/pe/netcdf/4.9.0.9/nvidia/23.3/lib" \
-DNETCDF_CXX_SHARED_LIB_DIR="/opt/cray/pe/netcdf/4.9.0.9/nvidia/23.3/lib" \
-DLIBCONFIG_C_HEADERS_DIR="/global/homes/t/tyounkin/code/libconfigBuild/gnu/include" \
-DLIBCONFIG_CXX_HEADERS_DIR="/global/homes/t/tyounkin/code/libconfigBuild/gnu/include" \
-DLIBCONFIG_C_LIB_DIR="/global/homes/t/tyounkin/code/libconfigBuild/gnu/lib" \
-DLIBCONFIG_CXX_LIB_DIR="/global/homes/t/tyounkin/code/libconfigBuild/gnu/lib" \
&> $2/cmake_output.txt ;

cmake --build $2 -- -j &> $2/build_output.txt
