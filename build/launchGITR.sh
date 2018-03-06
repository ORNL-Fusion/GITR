#!/bin/sh
echo "Running shell script for GITR mpirun"
#mpirun -n 1 -ppn 24 ${HOME}/gitr/build/GITR
srun -n 1 -c 24 $GITR_PATH/build/GITR
exit 0
