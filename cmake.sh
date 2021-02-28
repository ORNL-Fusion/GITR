rm -rf CMakeFiles
rm CMakeCache.txt
rm CTestTestfile.cmake
rm Makefile
rm cmake_install.cmake
rm install_manifest.txt
rm *test*
rm GITR
rm *.a
# Script to build GITR
cmake -DGITR_USE_CUDA=0 \
      -DCUDA_TOOLKIT_ROOT_DIR=/home/dg6/spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.4.0/cuda-10.2.89-hzby2dvk7t7hh4tbqtx66uexreoc2ih3 \
.. && make -j 4
