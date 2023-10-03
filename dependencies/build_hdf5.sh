wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.gz
tar -xzvf hdf5-1.12.2.tar.gz
cd hdf5-1.12.2/
mkdir ../hdf5build
./configure --with-zlib=$GITR_TOP_LEVEL/zlibbuild --prefix=$GITR_TOP_LEVEL/hdf5build
make -j
#make check
make install
make check-install
