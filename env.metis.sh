#!/bin/bash
module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load cray-netcdf
module load boost

export LD_LIBRARY_PATH=/home/tyounkin/code/libconfig/gnu/lib:$LD_LIBRARY_PATH
