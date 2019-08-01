#ifndef _SURFACES_
#define _SURFACES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
//#include <cmath>
#include "math.h"
#include <stdio.h>
//#include <vector>
#include "array.h"
//#include "managed_allocation.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#include <random>

//CUDA_CALLABLE_MEMBER

class Surfaces : public ManagedAllocation {
public: 
  int nSurfaces;  
  int nE;
  int nA;
  float E0;
  float E;
  float A0;
  float A;
  float dE;
  float dA;
  sim::Array<int> sumParticlesStrike;
  sim::Array<float> gridE;
  sim::Array<float> gridA;
  sim::Array<float> sumWeightStrike;
  sim::Array<float> grossDeposition;
  sim::Array<float> grossErosion;
  sim::Array<float> aveSputtYld;
  sim::Array<int> sputtYldCount;
  sim::Array<float> energyDistribution;
  sim::Array<float> sputtDistribution;
  sim::Array<float> reflDistribution;

  CUDA_CALLABLE_MEMBER
  
  void setSurface(int nE, float E0, float E, int nA, float A0, float A) {
    this->nE = nE;
    this->E0 = E0;
    this->E = E;
    this->nA = nA;
    this->A0 = A0;
    this->A = A;
    this->dE = (E - E0) / static_cast<float>(nE);
    this->dA = (A - A0) / static_cast<float>(nA);
    for (int i = 0; i < nE; i++) {
      this->gridE[i] = E0 + static_cast<float>(i) * dE;
    }

    for (int i = 0; i < nA; i++) {
      this->gridA[i] = A0 + static_cast<float>(i) * dA;
    }
  };

  CUDA_CALLABLE_MEMBER
  Surfaces(std::size_t nS,std::size_t nE, std::size_t nA) :
   energyDistribution{nS*nE*nA,0.0},sputtDistribution{nS*nE*nA,0.0},
   reflDistribution{nS*nE*nA,0.0}, gridE{nE,0.0}, gridA{nA,0.0},
 sumWeightStrike{nS,0.0}, sumParticlesStrike{nS,0}, grossDeposition{nS,0.0},
  grossErosion{nS,0.0}, aveSputtYld{nS,0.0}, sputtYldCount{nS,0} {};   

};

#endif
