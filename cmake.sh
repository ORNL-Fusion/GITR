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
#      -DCMAKE_CXX_COMPILER=/usr/lib/llvm-12/bin/clang \
#      -DCMAKE_C_COMPILER=/usr/lib/llvm-12/bin/clang \
cmake -DGITR_USE_CUDA=0 \
      -DMPIEXEC_EXECUTABLE=/home/dg6/spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.4.0/openmpi-3.1.5-s4r4ogvq4koshpy62bjvslz7bu33upxu/bin/mpiexec \
      -DGITR_USE_MPI=0 \
      -DMPI_CXX=/home/dg6/spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.4.0/openmpi-3.1.5-s4r4ogvq4koshpy62bjvslz7bu33upxu/lib/libmpi_cxx.so \
      -DCMAKE_CXX_COMPILER=clang-12\
      -DCMAKE_C_COMPILER=clang-12\
      -DCUDA_TOOLKIT_ROOT_DIR=/home/dg6/spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.4.0/cuda-10.2.89-hzby2dvk7t7hh4tbqtx66uexreoc2ih3 \
.. && make -j 4
