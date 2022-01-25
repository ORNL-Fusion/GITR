cmake -S /home/5n4/GITR -B /home/5n4/build_gitr \
&> /home/5n4/build_gitr/build_log.txt ;

cmake --build /home/5n4/build_gitr -- -j \
&>> /home/5n4/build_gitr/build_log.txt
