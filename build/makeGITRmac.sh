#!/bin/bash
source ../env.mac.sh
#CMAKE_C_COMPILER="/usr/local/Cellar/llvm/5.0.1/bin/clang"

#sCMAKE_CXX_COMPILER="/usr/local/Cellar/llvm/5.0.1/bin/clang++"
#OPENMP_LIBRARIES="/usr/local/Cellar/llvm/5.0.1/lib"
#OPENMP_INCLUDES="/usr/local/Cellar/llvm/5.0.1/include"
    #-DOpenMP_C=${CMAKE_C_COMPILER} \
    #-DOpenMP_CXX=${CMAKE_CXX_COMPILER} \
    #-DOpenMP_CXX_LIB_NAMES=libomp \
    #-DOpenMP_C_LIB_NAMES=libomp \
    #-DOpenMP_C_FLAGS="-Wno-unused-command-line-argument -fopenmp" \
    #-DOpenMP_CXX_FLAGS=" -Wno-unused-command-line-argument -fopenmp" \
    #-DOpenMP_libomp_LIBRARY=/usr/local/lib/libomp.dylib \
    #-DOpenMP_libgomp_LIBRARY=/usr/local/lib/libgomp \
    #-DOpenMP_libiomp5_LIBRARY=/usr/local/lib/libiomp5 \
    #    -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang  \
    #    -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ \
    
/opt/local/bin/cmake -DTHRUST_INCLUDE_DIR=$THRUST_DIR \
		     -DCMAKE_CXX_FLAGS="-g -Wall -Wextra -Wpedantic -Wno-unused-parameter -Wno-float-conversion -Wno-unused-variable -Wno-unused-but-const-variable" \
		     -DCMAKE_C_FLAGS="-g -Wall -Wextra -Wpedantic -Wno-unused-parameter -Wno-float-conversion -Wno-unused-variable -Wno-unused-but-const-variable" \
		     -DOpenMP_C_FLAGS="-fopenmp=libomp" \
		     -DOpenMP_C_LIB_NAMES="libomp" \
		     -DOpenMP_libomp_LIBRARY="/opt/local/lib/libomp/libomp.dylib" \
		     -DOpenMP_CXX_FLAGS="-fopenmp=libomp" \
		     -DOpenMP_CXX_LIB_NAMES="libomp" \
    -DNETCDF_CXX_INCLUDE_DIR=$NETCDFCXX4INCLUDE \
    -DNETCDF_CXX_LIBRARY=$NETCDFLIB_CPP \
    -DNETCDF_DIR=$NETCDFDIR \
    -DNETCDF_INCLUDE_DIR=$NETCDFINCLUDE \
    -DNETCDF_LIBRARY=$NETCDFLIB \
    -DNETCDF_CXX_INCLUDE_DIR=$NETCDFCXX4INCLUDE \
    -DLIBCONFIGPP_INCLUDE_DIR=$LIBCONFIGPP_INCLUDE \
    -DUSE_CUDA=0 \
    -DUSE_MPI=1 \
    -DUSEIONIZATION=0 \
    -DUSERECOMBINATION=0 \
    -DUSEPERPDIFFUSION=0 \
    -DUSECOULOMBCOLLISIONS=1 \
    -DUSEFRICTION=1 \
    -DUSEANGLESCATTERING=1 \
    -DUSEHEATING=1 \
    -DUSETHERMALFORCE=0 \
    -DUSESURFACEMODEL=0 \
    -DUSESHEATHEFIELD=0 \
    -DBIASED_SURFACE=0 \
    -DUSEPRESHEATHEFIELD=0 \
    -DBFIELD_INTERP=2 \
    -DLC_INTERP=0 \
    -DGENERATE_LC=0 \
    -DEFIELD_INTERP=0 \
    -DPRESHEATH_INTERP=0 \
    -DDENSITY_INTERP=2 \
    -DTEMP_INTERP=2 \
    -DFLOWV_INTERP=0 \
    -DGRADT_INTERP=0 \
    -DODEINT=0 \
    -DFIXEDSEEDS=1 \
    -DPARTICLESEEDS=1 \
    -DGEOM_TRACE=1 \
    -DGEOM_HASH=0 \
    -DGEOM_HASH_SHEATH=0 \
    -DPARTICLE_TRACKS=1 \
    -DPARTICLE_SOURCE_SPACE=0 \
    -DPARTICLE_SOURCE_ENERGY=0 \
    -DPARTICLE_SOURCE_ANGLE=0 \
    -DPARTICLE_SOURCE_FILE=1 \
    -DSPECTROSCOPY=2 \
    -DUSE3DTETGEOM=0 \
    -DUSECYLSYMM=1 \
    -DUSEFIELDALIGNEDVALUES=0 \
    -DFLUX_EA=1 \
    -DFORCE_EVAL=0 \
    -DUSE_SORT=0 \
    -DCHECK_COMPATIBILITY=1 \
    ..
    #-DCMAKE_CXX_FLAGS="-g -Wall -Wextra -Wpedantic -Wconversion -Wno-unused-parameter -Wno-float-conversion -Wno-unused-variable -Wno-unused-but-set-variable" \
