#!/bin/bash
module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
#module load cray-netcdf
module load cray-parallel-netcdf
module load boost/1.63
module load cray-mpich
#module load libconfig/gcc/64/1.5 
#module load git
#module load mpich/ge/gcc/64/3.1
##module load cuda75/toolkit/7.5.18
#module load cuda70/toolkit/7.0.28 
#module load slurm
##module load netcdf4
#PATH=$PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin
##prepend-path     PKG_CONFIG_PATH /home/dg6/code/netcdf/gnu/lib/pkgconfig 
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin
#
export LIBCONFIGDIR=/global/homes/t/tyounkin/code/libconfigBuild/gnu
export LIBCONFIGLIB=$LIBCONFIGDIR/lib
export LIBCONFIGPP_LIBRARIES=lconfig++
export LIBCONFIGPP_LIBRARY=lconfig++
export LIBCONFIGPP_STATIC_LIBRARY=
export LIBCONFIG_INCLUDE_DIR=$LIBCONFIGDIR/include
export LIBCONFIGPP_INCLUDE_DIR=$LIBCONFIGDIR/include
export LIBCONFIG_LIBRARY=lconfig
#LD_LIBRARY_PATH=LD_RUN_PATH
export THRUST_INCLUDE_DIRS=/global/homes/t/tyounkin/code/thrust/include
export NETCDF_CXX_INCLUDE_DIR=$NETCDF_DIR/include
export NETCDF_CXX_LIBRARY=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIR=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIRS=$NETCDF_DIR/include
export NETCDF_LIBRARIES=$NETCDF_DIR/lib
export NETCDF_LIBRARY=$NETCDF_DIR/lib
#export NETCDF=$NETCDF_DIR
#export NETCDF_CXX_INCLUDE_DIR=$NETCDF/include
#export NETCDF_CXX_LIBRARY=
#export NETCDF_INCLUDE_DIR=$NETCDF/include
#export NETCDF_LIBRARY
#export  NETCDFCXX4=/home/tqd/code/netcdfcxx  
#export  NETCDFDIR=/home/tqd/code/netcdfBuild/lib 
#export  NETCDFCXX4DIR=/home/tqd/code/netcdfcxx/lib 
#export  NETCDFINCLUDE=/home/tqd/code/netcdfBuild/include 
#export  NETCDFCXX4INCLUDE=/home/tqd/code/netcdfcxx/include 
#export  NETCDFLIB=netcdf 
#export  NETCDFLIB_CPP=netcdf_c++
#export NETCDF_LIBRARIES=$NETCDF/lib
#export NETCDF_INCLUDE_DIRS=$NETCDF/include
#
##boost vars
#export Boost_INCLUDE_DIRS=~/code/boost/include
#export Boost_LIBRARY_DIRS=/global/homes/t/tyounkin/code/boostBuild/lib
#export BOOST_ROOT=~/code/boost
#export BOOST_INCLUDEDIR=~/code/boost/include
#export BOOST_LIBRARYDIR=~/code/boost/lib
#
#export PYTHONPATH=/home/tqd/code/netcdf4-python/build
#export MPI_C_LIBRARIES=mpich
#export MPI_C_INCLUDE_PATH=/opt/cray/pe/mpt/7.6.2/gni/mpich-gnu/5.1/include
#export MPI_CXX_LIBRARIES=mpichcxx
#export MPI_CXX_INCLUDE_PATH=/opt/cray/pe/mpt/7.6.2/gni/mpich-gnu/5.1/include
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBCONFIGLIB:/opt/cray/pe/mpt/7.6.2/gni/mpich-gnu/5.1/lib
#export LD_LIBRARY_PATH

#export OMPI_ROOT=/project/projectdirs/atom/users/elwasif/ompi/install_4.0
#export PATH=$OMPI_ROOT/bin:$PATH
#export LD_LIBRARY_PATH=$OMPI_ROOT/lib64:$OMPI_ROOT/lib:$LD_LIBRARY_PATH
#export MANPATH=$OMPI_ROOT/share/man:$MANPATH
export ATOM=/project/projectdirs/atom

################
# Edison specfic
################

export ATOM_EDISON=$ATOM/atom-install-edison
export GITR_PATH=$ATOM/atom-install-edison/GITR
export FTRIDYN_PATH=$ATOM/atom-install-cori/fractal-tridyn
export PYTHONPATH=$GITR_PATH/ftridyn:/global/homes/t/tyounkin/code/mpi4pyBuild/lib.linux-x86_64-2.7:$PYTHONPATH:$FTRIDYN_PATH/utils:$GITR_PATH/python:/project/projectdirs/atom/users/tyounkin/libconfPython:/project/projectdirs/atom/users/tyounkin/libconfPython/lib/python2.7/site-packages/:/global/homes/t/tyounkin/code/netcdfPython/lib.linux-x86_64-2.7/
