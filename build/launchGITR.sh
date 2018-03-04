#!/bin/sh
echo "Running shell script for GITR mpirun"
mpirun -n 8 -ppn 4 ${HOME}/gitr/build/GITR
exit 0
