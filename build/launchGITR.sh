#!/bin/sh
echo "Running shell script for GITR mpirun"
#mpirun -n 2 /project/projectdirs/atom/atom-install-edison/GITR/build/GITR
srun -n 10 -c 24 $GITR_PATH/build/GITR
exit 0
