#!/bin/sh
echo "Running shell script for GITR mpirun"
mpirun -n 8 -ppn 4 ${HOME}/gitr/build/GITR
#mpirun -n 2 /project/projectdirs/atom/atom-install-edison/GITR/build/GITR
#srun -n 1 -c 24 $GITR_PATH/build/GITR
exit 0
