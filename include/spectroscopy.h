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

#if USE_OPENMP == 1
#include "omp.h"
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
#if USE_CUDA >0
__device__ double atomicAdd1(gitr_precision* address, gitr_precision val);
#endif

struct spec_bin { 
    Flags *flags;
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
    int cylsymm;
    int spectroscopy;

    spec_bin(Flags* _flags, Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, gitr_precision *_gridX,gitr_precision *_gridY,gitr_precision *_gridZ,
           gitr_precision* _bins, gitr_precision _dt, int cylsymm_, int spectroscopy_ );

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const;

};

#endif
