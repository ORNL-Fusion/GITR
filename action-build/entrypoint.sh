#!/bin/sh -l

sh -c "echo directory $PWD"
sh -c "ls /home/runner/work/GITR-1/GITR-1"
sh -c "cd build"
sh -c "echo now directory $PWD"
sh -c "ls /github/workspace"

sh -c "build/makeGITR.sh"
