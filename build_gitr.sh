cmake -S /home/5n4/GITR -B /home/5n4/build \
-DGITR_USE_PARTICLE_SOURCE_FILE=0 \
-DGITR_USE_SURFACE_POTENTIAL=0 \
-DGITR_USE_3DTET_GEOM=0 \
-DGITR_USE_CUDA=0 -DGITR_USE_MPI=1 \
-DCMAKE_BUILD_TYPE=Debug \
-DGITR_SPECTROSCOPY=2 \
-DGITR_USE_CYLSYMM=0 \
-DGITR_BFIELD_INTERP=0 \
&> /home/5n4/build/build_log.txt;

cmake --build /home/5n4/build -- -j &>> /home/5n4/build/build_log.txt
