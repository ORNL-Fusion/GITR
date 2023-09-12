#ifndef _SURFACES_
#define _SURFACES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include <cstdlib>
#include <string>

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#include <random>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

//CUDA_CALLABLE_MEMBER
class Surfaces : public ManagedAllocation {
public: 
  int nSurfaces;  
  int nSpecies;
  int nE;
  int nA;
  gitr_precision E0;
  gitr_precision E;
  gitr_precision A0;
  gitr_precision A;
  gitr_precision dE;
  gitr_precision dA;
  sim::Array<gitr_precision> sumParticlesStrike;
  sim::Array<gitr_precision> gridE;
  sim::Array<gitr_precision> gridA;
  sim::Array<gitr_precision> sumWeightStrike;
  sim::Array<gitr_precision> grossDeposition;
  sim::Array<gitr_precision> grossErosion;
  sim::Array<gitr_precision> aveSputtYld;
  sim::Array<gitr_precision> sputtYldCount;
  sim::Array<gitr_precision> energyDistribution;
  sim::Array<gitr_precision> sputtDistribution;
  sim::Array<gitr_precision> reflDistribution;

  CUDA_CALLABLE_MEMBER
  void setSurface(int nE, gitr_precision E0, gitr_precision E, int nA, gitr_precision A0, gitr_precision A) {
    this->nE = nE;
    this->E0 = E0;
    this->E = E;
    this->nA = nA;
    this->A0 = A0;
    this->A = A;
    this->dE = (E - E0) / static_cast<gitr_precision>(nE);
    this->dA = (A - A0) / static_cast<gitr_precision>(nA);
    for (int i = 0; i < nE; i++) {
      this->gridE[i] = E0 + static_cast<gitr_precision>(i) * dE;
    }
    for (int i = 0; i < nA; i++) {
      this->gridA[i] = A0 + static_cast<gitr_precision>(i) * dA;
    }
  };

CUDA_CALLABLE_MEMBER
Surfaces(std::size_t nS, std::size_t nSpecies,
         std::size_t nE, std::size_t nA) :
    sumParticlesStrike{nS * nSpecies, 0.0},
    gridE{nE,0.0},
    gridA{nA,0.0},
    sumWeightStrike{nS * nSpecies,0.0},
    grossDeposition{nS * nSpecies,0.0},
    grossErosion{nS * nSpecies,0.0},
    aveSputtYld{nS * nSpecies,0.0},
    sputtYldCount{nS * nSpecies,0},
    energyDistribution{nS * nE * nA, 0.0},
    sputtDistribution{nS * nE * nA,0.0},
    reflDistribution{nS * nE * nA,0.0} {};   
};

#endif // _SURFACES_