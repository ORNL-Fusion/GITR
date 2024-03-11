export GITR_DEPENDENCIES=~/GITR/dependencies/
cmake -S $1 -B $2 \
-DNETCDF_C_HEADERS_DIR="~/GITR/dependencies/netcdfcbuild/include" \
-DNETCDF_CXX_HEADERS_DIR="~/GITR/dependencies/netcdfcxxbuild/include" \
-DNETCDF_C_SHARED_LIB_DIR="~/GITR/dependencies/netcdfcbuild/lib" \
-DNETCDF_CXX_SHARED_LIB_DIR="~/GITR/dependencies/netcdfcxxbuild/lib" \
-DLIBCONFIG_C_HEADERS_DIR="~/GITR/dependencies/libconfigbuild/include" \
-DLIBCONFIG_CXX_HEADERS_DIR="~/GITR/dependencies/libconfigbuild/include" \
-DLIBCONFIG_C_LIB_DIR="~/GITR/dependencies/libconfigbuild/lib" \
-DLIBCONFIG_CXX_LIB_DIR="~/GITR/dependencies/libconfigbuild/lib" \
.
