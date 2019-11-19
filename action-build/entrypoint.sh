#!/bin/sh -l

sh -c "echo directory $PWD"
sh -c "ls /home/runner/work/GITR-1/GITR-1"
cd GITR-1/build
sh -c "echo now directory $PWD"
sh -c "ls /github/workspace"

./makeGITR.sh
