#ifndef _HASHGEOM_
#define _HASHGEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include "Boundary.h"
#include "boris.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#endif

#ifdef __GNUC__ 
#include <random>
#endif

#include "interpRateCoeff.hpp"
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct hashGeom {
   //int k;
   int nLines; 
   int nHashes;
   Boundary* boundary;
   gitr_precision* x;
   gitr_precision* y;
   gitr_precision* z;
   int* n_closeGeomElements;
   //gitr_precision* minDist;
   int* closeGeom;
   int* nR;
   int* nY;
   int* nZ;
   int use_3d_geom;


   hashGeom( int _nLines,int _nHashes,
                Boundary* _boundary,
                gitr_precision* _x,
                gitr_precision* _y, 
                gitr_precision* _z, 
                int* _n_closeGeomElements,//gitr_precision *_minDist,
                int *_closeGeom,
                int* _nR, int* _nY, int* _nZ, int use_3d_geom_ );
    
   CUDA_CALLABLE_MEMBER_DEVICE 
   void operator()(std::size_t indx); 
};

#endif
