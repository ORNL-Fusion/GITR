rm -rf CMakeFiles
rm CMakeCache.txt
rm CTestTestfile.cmake
rm Makefile
rm cmake_install.cmake
rm install_manifest.txt
rm *test*
rm GITR
rm *.a

cmake -DGITR_USE_CUDA=0 \
.. && make -j 4
