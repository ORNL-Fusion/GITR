module load cmake
module load cuda/10.0
module load openmpi/3.0.0/1
export CODE_PATH=$HOME/barn/code
#$PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/tqd/code/netcdfBuild/bin:/home/tqd/code/netcdfcxx/bin

export  NETCDF=$CODE_PATH/netcdfCBuild 
export  NETCDFCXX4=$CODE_PATH/netcdfCXXBuild  
export  NETCDFDIR=$CODE_PATH/netcdfCBuild/lib 
export  NETCDFCXX4DIR=$CODE_PATH/netcdfCXXBuild/lib 
export  NETCDFINCLUDE=$CODE_PATH/netcdfCBuild/include 
export  NETCDFCXX4INCLUDE=$CODE_PATH/netcdfCXXBuild/include 
export  NETCDFLIB=netcdf 
export  NETCDFLIB_CPP=netcdf_c++4
export NETCDF_LIBRARIES=$NETCDFCXX4DIR/$NETCDFLIB_CPP 
export NETCDF_INCLUDE_DIRS=$NETCDFCXX4INCLUDE
#boost vars
export Boost_INCLUDE_DIRS=$CODE_PATH/boostBuild/gnu/include
export Boost_LIBRARY_DIRS=$CODE_PATH/boostBuild/gnu/lib
export Boost_LIBRARIES_LIST=${$Boost_LIBRARY_DIRS/libboost_chrono.so}
export BOOST_ROOT=$CODE_PATH/boostBuild/gnu/boost
export BOOST_INCLUDEDIR=$CODE_PATH/boostBuild/gnu/include
export BOOST_LIBRARYDIR=$CODE_PATH/boostBuild/gnu/lib
export MPI_C_LIBRARIES=/usr/mpi/gcc/openmpi-1.8.4/lib64
export LIBCONFIGPPINCLUDE=$CODE_PATH/libconfigBuild/include
export LIBCONFIGPPLIB=$CODE_PATH/libconfigBuild/lib
PYTHONPATH=$PYTHONPATH:/home/tqd/gitr/python/
export PYTHONPATH=$PYTHONPATH:/home/tqd/code/netcdf4-python/build
