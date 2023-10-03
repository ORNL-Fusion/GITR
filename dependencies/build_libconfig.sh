wget https://github.com/hyperrealm/libconfig/archive/v1.7.2.tar.gz
tar -xvzf v1.7.2.tar.gz
cd libconfig-1.7.2/
mkdir ../libconfigbuild
autoreconf
./configure --prefix=$GITR_TOP_LEVEL/libconfigbuild
make
make check
make install
