wget https://ftp.gnu.org/gnu/m4/m4-1.4.19.tar.gz
tar -xzvf m4-1.4.19.tar.gz
cd m4-1.4.19/
mkdir ../m4build
./configure --prefix=$GITR_TOP_LEVEL/m4build
make
make install
