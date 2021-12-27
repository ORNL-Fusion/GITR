cmake -S /home/5n4/GITR -B /home/5n4/build_gitr -G Ninja \
-DGITR_USE_CUDA=0 \
-DGITR_USE_PRE_SHEATH_EFIELD=1 \
-DGITR_USE_SHEATH_EFIELD=0 \
&> /home/5n4/build_gitr/build_log.txt ;

#cmake --build /home/5n4/build_gitr -- -j 0 &>> /home/5n4/build_gitr/build_log.txt
