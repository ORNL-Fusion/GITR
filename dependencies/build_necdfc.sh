wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.7.4.tar.gz
tar -xzvf v4.7.4.tar.gz
cd netcdf-c-4.7.4/
mkdir ../netcdfcbuild
export CPPFLAGS="-I$GITR_TOP_LEVEL/zlibbuild/include -I$GITR_TOP_LEVEL/hdf5build/include"
export LDFLAGS="-L$GITR_TOP_LEVEL/zlibbuild/lib -L$GITR_TOP_LEVEL/hdf5build/lib"
./configure --disable-dap --prefix=$GITR_TOP_LEVEL/netcdfcbuild
make
#make check
make install
