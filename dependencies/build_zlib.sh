wget https://zlib.net/current/zlib.tar.gz
tar -xzvf zlib.tar.gz
cd zlib-1.3.1
mkdir ../zlibbuild
./configure --prefix=$GITR_TOP_LEVEL/zlibbuild
make
make install
