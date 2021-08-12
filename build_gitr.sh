/home/5n4/cmake-3.21.1_build/bin/cmake -S /home/5n4/GITR -B /home/5n4/build \
-DGITR_USE_CUDA=1 \
&> /home/5n4/build/build_log.txt ;
/home/5n4/cmake-3.21.1_build/bin/cmake --build /home/5n4/build -- -j \
&>> /home/5n4/build/build_log.txt
