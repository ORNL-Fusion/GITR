#!/bin/sh
echo "Running shell script for Fractal Tridyn Input file FTridyn.IN"
$FTRIDYN_PATH/src/a.out < $1
exit 0
