#!/bin/bash
module load libconfig/gcc/64/1.5 
module load git
#module load cuda75/toolkit/7.5.18
module load cuda70/toolkit/7.0.28 
module load slurm
#module load netcdf4
PATH=$PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin
#prepend-path     PKG_CONFIG_PATH /home/dg6/code/netcdf/gnu/lib/pkgconfig 
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin

export  NETCDF=/home/tqd/code/netcdfBuild 
export  NETCDFCXX4=/home/tqd/code/netcdfcxx  
export  NETCDFDIR=/home/tqd/code/netcdfBuild/lib 
export  NETCDFCXX4DIR=/home/tqd/code/netcdfcxx/lib 
export  NETCDFINCLUDE=/home/tqd/code/netcdfBuild/include 
export  NETCDFCXX4INCLUDE=/home/tqd/code/netcdfcxx/include 
export  NETCDFLIB=netcdf 
export  NETCDFLIB_CPP=netcdf_c++ 
