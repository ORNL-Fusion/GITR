/home/5n4/cmake-3.21.1_build/bin/cmake -S /home/5n4/GITR -B /home/5n4/build_gitr \
-DGITR_USE_CUDA=1 \
-DGITR_USE_PRE_SHEATH_EFIELD=1 \
-DGITR_USE_SHEATH_EFIELD=0 \
&> /home/5n4/build_gitr/build_log.txt ;
/home/5n4/cmake-3.21.1_build/bin/cmake --build /home/5n4/build_gitr -- -j \
&>> /home/5n4/build_gitr/build_log.txt