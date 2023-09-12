#ifndef _SURFACE_
#define _SURFACE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include "Surfaces.h"
#include <cmath>
#include "boris.h"
#include "spectroscopy.h"


#if USE_OPENMP == 1
#include "omp.h"
#endif

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif


struct reflection {
    Particles * particles;
    const double dt;
    int nLines;
    Boundary * boundaryVector;
    Surfaces * surfaces;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    int flux_ea;
    int use_3d_geom;
    int cylsymm;
    int nspecies;

    reflection(Particles* _particles, double _dt,
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nLines,Boundary * _boundaryVector,
            Surfaces * _surfaces,
    int flux_ea_,
    int use_3d_geom_,
    int cylsymm_ ,
        int nspecies_);

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const;
    
};
#endif
