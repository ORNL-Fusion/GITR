#!/bin/sh -l

sh -c "echo directory $PWD"

cd GITR-1/build
sh -c "echo now directory $PWD"

./makeGITR.sh
