cmake -S $1 -B $2 -G Ninja &> $2/cmake_output.txt ;

cmake --build $2 -- -j 0 &> $2/build_output.txt
