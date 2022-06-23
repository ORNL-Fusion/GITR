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
    int nE_sputtRefCoeff;
    int nA_sputtRefCoeff;
    gitr_precision* A_sputtRefCoeff;
    gitr_precision* Elog_sputtRefCoeff;
    gitr_precision* spyl_surfaceModel;
    gitr_precision* rfyl_surfaceModel;
    int nE_sputtRefDistOut; 
    int nE_sputtRefDistOutRef; 
    int nA_sputtRefDistOut;
    int nE_sputtRefDistIn;
    int nA_sputtRefDistIn;
    gitr_precision* E_sputtRefDistIn;
    gitr_precision* A_sputtRefDistIn;
    gitr_precision* E_sputtRefDistOut;
    gitr_precision* E_sputtRefDistOutRef;
    gitr_precision* A_sputtRefDistOut;
    gitr_precision* energyDistGrid01;
    gitr_precision* energyDistGrid01Ref;
    gitr_precision* angleDistGrid01;
    gitr_precision* EDist_CDF_Y_regrid;
    gitr_precision* ADist_CDF_Y_regrid;
    gitr_precision* EDist_CDF_R_regrid;
    gitr_precision* ADist_CDF_R_regrid;
    int nEdist;
    gitr_precision E0dist;
    gitr_precision Edist;
    int nAdist;
    gitr_precision A0dist;
    gitr_precision Adist;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    int flux_ea;
    int use_3d_geom;
    int cylsymm;

    reflection(Particles* _particles, double _dt,
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nLines,Boundary * _boundaryVector,
            Surfaces * _surfaces,
    int _nE_sputtRefCoeff,
    int _nA_sputtRefCoeff,
    gitr_precision* _A_sputtRefCoeff,
    gitr_precision* _Elog_sputtRefCoeff,
    gitr_precision* _spyl_surfaceModel,
    gitr_precision* _rfyl_surfaceModel,
    int _nE_sputtRefDistOut,
    int _nE_sputtRefDistOutRef,
    int _nA_sputtRefDistOut,
    int _nE_sputtRefDistIn,
    int _nA_sputtRefDistIn,
    gitr_precision* _E_sputtRefDistIn,
    gitr_precision* _A_sputtRefDistIn,
    gitr_precision* _E_sputtRefDistOut,
    gitr_precision* _E_sputtRefDistOutRef,
    gitr_precision* _A_sputtRefDistOut,
    gitr_precision* _energyDistGrid01,
    gitr_precision* _energyDistGrid01Ref,
    gitr_precision* _angleDistGrid01,
    gitr_precision* _EDist_CDF_Y_regrid,
    gitr_precision* _ADist_CDF_Y_regrid, 
    gitr_precision* _EDist_CDF_R_regrid,
    gitr_precision* _ADist_CDF_R_regrid,
    int _nEdist,
    gitr_precision _E0dist,
    gitr_precision _Edist,
    int _nAdist,
    gitr_precision _A0dist,
    gitr_precision _Adist,
    int flux_ea_,
    int use_3d_geom_,
    int cylsymm_ );

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const;
    
};
#endif
