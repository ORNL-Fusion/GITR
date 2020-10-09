module swap PrgEnv-intel PrgEnv-gnu
module load cmake3/3.6.1   #3.2.3
#module load cudatoolkit
module load cray-hdf5
module load cray-netcdf
#module load boost/1.62.0
#module load libconfig


#export NETCDF_DIR=/opt/cray/netcdf/4.4.1.1.3/gnu/4.9/
export NETCDF_CXX_INCLUDE_DIR=$NETCDF_DIR/include
export NETCDF_CXX_LIBRARY=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIR=$NETCDF_DIR/lib
export NETCDF_INCLUDE_DIRS=$NETCDF_DIR/include
export NETCDF_LIBRARIES=$NETCDF_DIR/lib
export NETCDF_LIBRARY=$NETCDF_DIR/lib
export LIBCONFIGDIR=$MEMBERWORK/fus049/code/libconfigBuild
#export LIBCONFIGDIR=/ccs/home/tyounkin/code/libconfigBuild
export LIBCONFIGLIB=$LIBCONFIGDIR/lib
export LIBCONFIGPP_LIBRARIES=lconfig++

#boost vars
#export BOOST_ROOT=/ccs/home/tyounkin/code/boostBuild
export BOOST_ROOT=$MEMBERWORK/fus049/code/boostBuild
export Boost_INCLUDE_DIRS=$BOOST_ROOT/include/boost
export Boost_LIBRARY_DIRS=$BOOST_ROOT/lib
export BOOST_INCLUDEDIR=$BOOST_ROOT/include/boost
export BOOST_LIBRARYDIR=$BOOST_ROOT/lib

#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/atlas/scratch/tyounkin/fus049/code/boostBuild/lib
export THRUST_INCLUDE_DIRS=/ccs/home/tyounkin/code/thrust
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$THRUST_INCLUDE_DIRS:$LIBCONFIGLIB:$BOOST_LIBRARYDIR
