wget http://zlib.net/zlib-1.2.13.tar.gz
tar -xzvf zlib-1.2.13.tar.gz
cd zlib-1.2.13
mkdir ../zlibbuild
./configure --prefix=$GITR_TOP_LEVEL/zlibbuild
make
make install
