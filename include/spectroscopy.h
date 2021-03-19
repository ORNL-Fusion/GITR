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
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
#if USE_CUDA >0
__device__ double atomicAdd1(double* address, double val);
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
    double *bins;
    gitr_precision dt;

    spec_bin(Flags* _flags, Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, gitr_precision *_gridX,gitr_precision *_gridY,gitr_precision *_gridZ,
           double * _bins, gitr_precision _dt) : 
        flags(_flags), particlesPointer(_particlesPointer), nBins(_nBins),nX(_nX),nY(_nY), nZ(_nZ), gridX(_gridX),gridY(_gridY),gridZ(_gridZ), bins(_bins),
        dt(_dt) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const;

};

#endif
