#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Boundary.h"
#include "Particles.h"
#include "Surfaces.h"
#include "surfaceModel.h"
#include "boris.h"

#if USE_OPENMP == 1
#include "omp.h"
#endif

#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

__host__ __device__
gitr_precision
findT( gitr_precision x0, 
       gitr_precision x1,
       gitr_precision y0,
       gitr_precision y1,
       gitr_precision intersectionx );

  //template<int HOST=1>
struct geometry_check {
  Particles *particlesPointer;
  const int nLines;
  Boundary *boundaryVector;
  Surfaces *surfaces;
  gitr_precision dt;
  // int& tt;
  int nHashes;
  int *nR_closeGeom;
  int *nY_closeGeom;
  int *nZ_closeGeom;
  int *n_closeGeomElements;
  gitr_precision *closeGeomGridr;
  gitr_precision *closeGeomGridy;
  gitr_precision *closeGeomGridz;
  int *closeGeom;
  int nEdist;
  gitr_precision E0dist;
  gitr_precision Edist;
  int nAdist;
  gitr_precision A0dist;
  gitr_precision Adist;
  int flux_ea;
  int surface_model;
  int geom_hash;
  int use_3d_geom;
  int cylsymm;

  geometry_check(Particles *_particlesPointer, int _nLines,
                 Boundary *_boundaryVector, Surfaces *_surfaces, gitr_precision _dt,
                 int _nHashes, int *_nR_closeGeom, int *_nY_closeGeom,
                 int *_nZ_closeGeom, int *_n_closeGeomElements,
                 gitr_precision *_closeGeomGridr, gitr_precision *_closeGeomGridy,
                 gitr_precision *_closeGeomGridz, int *_closeGeom, int _nEdist,
                 gitr_precision _E0dist, gitr_precision _Edist, int _nAdist, gitr_precision _A0dist,
                 gitr_precision _Adist,
                 int flux_ea_,
                 int surface_model_,
                 int geom_hash_,
                 int use_3d_geom_,
                 int cylsymm_ );

  //CUDA_CALLABLE_MEMBER_DEVICE
 __host__  __device__
  void operator()(std::size_t indx) const;
};

#endif
