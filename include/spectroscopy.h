#ifndef _SPECTROSCOP_
#define _SPECTROSCOP_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <cmath>
#include <vector>
#include "config_interface.h"

#if USE_OPENMP == 1
#include "omp.h"
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

#if USE_CUDA >0
#include "atomic_add_1.h"
#endif

struct spec_bin { 
    class flags config_flags;
    Particles *particlesPointer;
    const int nBins;
    int nX;
    int nY;
    int nZ;
    gitr_precision *gridX;
    gitr_precision *gridY;
    gitr_precision *gridZ;
    gitr_precision *bins;
    gitr_precision dt;
    gitr_precision *bins_vx;
    gitr_precision *bins_vy;
    gitr_precision *bins_vz;
    gitr_precision *bins_E;

    spec_bin(class flags &f_init, Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, gitr_precision *_gridX,gitr_precision *_gridY,gitr_precision *_gridZ,
           gitr_precision* _bins, gitr_precision _dt,
           gitr_precision* _bins_vx,gitr_precision* _bins_vy,gitr_precision* _bins_vz, gitr_precision* _bins_E );

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const;

};

#endif
