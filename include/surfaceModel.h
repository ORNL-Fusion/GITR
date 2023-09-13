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

void processParticleHit(std::size_t indx) const;
void printMaterials(const std::string& incidentMaterial, const std::string& targetMaterial) const;
std::pair<std::string, std::string> getMaterials(int wallHit, std::size_t indx) const;

void reflect(Particles* particles, int indx, gitr_precision newWeight, gitr_precision eInterpVal, gitr_precision aInterpVal, 
             Boundary* boundaryVector, int wallHit, bool use_3d_geom, bool cylsymm,
             gitr_precision r10, gitr_precision* surfaceNormalVector) const;

void sputter(Boundary* boundaryVector, int wallHit,Particles* particles,int indx,gitr_precision aInterpVal,
    gitr_precision r10, gitr_precision newWeight,int nspecies,    bool use_3d_geom, bool cylsymm,
    const gitr_precision* surfaceNormalVector) const;

std::pair<gitr_precision, gitr_precision> computeIncidentParticleEnergyAngle(Particles* particles, int indx, int use_3d_geom, int cylsymm, gitr_precision* surfaceNormalVector) const;

gitr_precision getRandomValueForDevice(int indx) const;
bool sputteringEvent(gitr_precision randomReflector, gitr_precision sputtProb, gitr_precision totalYR) const;
bool reflectionEvent(gitr_precision randomReflector, gitr_precision sputtProb, gitr_precision totalYR) const;

};


#endif
