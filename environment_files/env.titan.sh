module swap PrgEnv-pgi/5.2.82 PrgEnv-gnu
module load cmake3
module load cudatoolkit
module load cray-netcdf
#module load boost/1.62.0
module load libconfig


#export NETCDF_DIR=/opt/cray/netcdf/4.4.1.1.3/gnu/4.9/
export NETCDF_CXX_INCLUDE_DIR=$NETCDF_DIR/include
export NETCDF_CXX_LIBRARY=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIR=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIRS=$NETCDF_DIR/include
export NETCDF_LIBRARIES=$NETCDF_DIR/lib
export NETCDF_LIBRARY=$NETCDF_DIR/lib
export LIBCONFIGDIR=/sw/xk6/libconfig/1.4.9/sles11.1_gnu4.3.4
export LIBCONFIGLIB=$LIBCONFIGDIR/lib
export LIBCONFIGPP_LIBRARIES=lconfig++

#boost vars
export BOOST_ROOT=/lustre/atlas/scratch/tyounkin/fus049/code/boostBuild
#export BOOST_ROOT=/ccs/home/tyounkin/code/boostBuild
export Boost_INCLUDE_DIRS=$BOOST_ROOT/include
export Boost_LIBRARY_DIRS=$BOOST_ROOT/lib
export BOOST_INCLUDEDIR=$BOOST_ROOT/include
export BOOST_LIBRARYDIR=$BOOST_ROOT/lib

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/atlas/scratch/tyounkin/fus049/code/boostBuild/lib
