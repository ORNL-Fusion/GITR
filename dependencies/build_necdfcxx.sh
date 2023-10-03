wget https://github.com/Unidata/netcdf-cxx4/archive/refs/tags/v4.3.1.tar.gz
tar -xzvf v4.3.1.tar.gz
cd netcdf-cxx4-4.3.1/
mkdir ../netcdfcxxbuild
export PATH=$PATH:/home/tqd/Code/netcdfcbuild/bin
export LDFLAGS="-L$GITR_TOP_LEVEL/zlibbuild/lib -L$GITR_TOP_LEVEL/hdf5build/lib -L$GITR_TOP_LEVEL/netcdfcbuild/lib"
export CPPFLAGS="-I$GITR_TOP_LEVEL/zlibbuild/include -I$GITR_TOP_LEVEL/hdf5build/include -I$GITR_TOP_LEVEL/netcdfcbuild/include"
./configure --prefix=$GITR_TOP_LEVEL/netcdfcxxbuild
make
#make check
make install
