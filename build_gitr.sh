cmake -S $1 -B $2 -G Ninja \
-DGITR_USE_CUDA=1 \
-DGITR_USE_PRE_SHEATH_EFIELD=1 \
-DGITR_USE_SHEATH_EFIELD=0 \
&> $2/cmake_output.txt ;

cmake --build $2 -- -j 0 &> $2/build_output.txt
