#!/bin/sh -l

sh -c "directory $PWD"

cd build
sh -c "now directory $PWD"

./makeGITR.sh
