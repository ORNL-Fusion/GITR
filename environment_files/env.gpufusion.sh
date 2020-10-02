#!/bin/bash
module load libconfig/gcc/64/1.5 
module load git
module load mpich
#module load openmpi
#module load cuda75/toolkit/7.5.18
module load cuda91
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

#boost vars
export Boost_INCLUDE_DIRS=~/code/boost/include
export Boost_LIBRARY_DIRS=~/code/boost/lib
export BOOST_ROOT=~/code/boost
export BOOST_INCLUDEDIR=~/code/boost/include
export BOOST_LIBRARYDIR=~/code/boost/lib
export MPI_C_LIBRARIES=/usr/mpi/gcc/openmpi-1.8.4/lib64

PYTHONPATH=$PYTHONPATH:/home/tqd/gitr/python/
export PYTHONPATH=$PYTHONPATH:/home/tqd/code/netcdf4-python/build
