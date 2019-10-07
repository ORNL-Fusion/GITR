#!/bin/bash
module swap PrgEnv-intel/6.0.5 PrgEnv-gnu/6.0.5
module load cray-netcdf
module load openmpi

export GITR_BASE_DIR=/project/projectdirs/m1709/psi-install-cori
export LIBCONFIGDIR=$GITR_BASE_DIR/libconfig
export LIBCONFIGLIB=$LIBCONFIGDIR/lib
export LIBCONFIGPP_LIBRARIES=lconfig++
export LIBCONFIGPP_LIBRARY=$LIBCONFIGLIB/libconfig++.so
export LIBCONFIGPP_STATIC_LIBRARY=
export LIBCONFIG_INCLUDE_DIR=$LIBCONFIGDIR/include
export LIBCONFIGPP_INCLUDE_DIR=$LIBCONFIGDIR/include
export LIBCONFIG_LIBRARY=lconfig
export THRUST_INCLUDE_DIRS=$GITR_BASE_DIR/thrust
export NETCDF_CXX_BASEDIR=/opt/cray/pe/netcdf/4.6.1.3/GNU/8.2
export NETCDF_CXX_INCLUDE_DIR=$NETCDF_CXX_BASEDIR/include
export NETCDF_CXX_LIBRARY=$NETCDF_CXX_BASEDIR/lib
export NETCDF_INCLUDE_DIR=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIRS=$NETCDF_DIR/include
export NETCDF_LIBRARIES=$NETCDF_DIR/lib
export NETCDF_LIBRARY=$NETCDF_DIR/lib

#export PATH=/global/homes/t/tyounkin/code/python2.7.14/bin/:$PATH
#export PYTHONPATH=/global/homes/t/tyounkin/code/libconfPython/lib/python2.7/site-packages/:$PYTHONPATH
#export PYTHONPATH=/global/homes/t/tyounkin/code/scipyBuild/lib/python2.7/site-packages:$PYTHONPATH
#export PYTHONPATH=/global/homes/t/tyounkin/atomIPS/atom-install-cori/GITR/python:$PYTHONPATH
#export PYTHONPATH=/global/homes/t/tyounkin/code/netcdfPython/lib.linux-x86_64-2.7:$PYTHONPATH
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBCONFIGLIB
#export LD_LIBRARY_PATH
