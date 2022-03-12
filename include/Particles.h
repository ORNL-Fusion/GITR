#ifndef _PARTICLES_
#define _PARTICLES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include <libconfig.h++>
#include "utils.h"

#include <cstdlib>
#include <stdio.h>

#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#endif

#include <random>
#include <math.h>
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

class Particles : public ManagedAllocation {
public:
  std::size_t nParticles;
  sim::Array<int> index;
  sim::Array<gitr_precision> x;
  sim::Array<gitr_precision> y;
  sim::Array<gitr_precision> z;
  sim::Array<gitr_precision> xprevious;
  sim::Array<gitr_precision> yprevious;
  sim::Array<gitr_precision> zprevious;
  sim::Array<gitr_precision> v;
  sim::Array<gitr_precision> vx;
  sim::Array<gitr_precision> vy;
  sim::Array<gitr_precision> vz;
  sim::Array<gitr_precision> Z;
  sim::Array<gitr_precision> amu;
  sim::Array<gitr_precision> charge;
  sim::Array<gitr_precision> newVelocity;
  sim::Array<gitr_precision> nu_s;
  sim::Array<gitr_precision> vD;
  sim::Array<int> tt;
  sim::Array<int> hasLeaked;
  sim::Array<gitr_precision> leakZ;
#ifdef __CUDACC__
  sim::Array<curandState> stream_ionize;
#else
  sim::Array<std::mt19937> stream_ionize;
#endif

  sim::Array<gitr_precision> hitWall;
  sim::Array<int> surfaceHit;
  sim::Array<int> firstCollision;
  sim::Array<gitr_precision> transitTime;
  sim::Array<gitr_precision> distTraveled;
  sim::Array<int> wallIndex;
  sim::Array<gitr_precision> perpDistanceToSurface;
  sim::Array<gitr_precision> test;
  sim::Array<gitr_precision> test0;
  sim::Array<gitr_precision> test1;
  sim::Array<gitr_precision> test2;
  sim::Array<gitr_precision> test3;
  sim::Array<gitr_precision> test4;
  sim::Array<gitr_precision> distanceTraveled;
  sim::Array<gitr_precision> weight;
  sim::Array<gitr_precision> firstIonizationZ;
  sim::Array<gitr_precision> firstIonizationT;
  sim::Array<gitr_precision> dt;
  sim::Array<gitr_precision> time;
  sim::Array<bool> advance;

  CUDA_CALLABLE_MEMBER
  void setParticle(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                   gitr_precision Ex, gitr_precision Ey, gitr_precision Ez,
                   gitr_precision Z, gitr_precision amu, gitr_precision charge);

  CUDA_CALLABLE_MEMBER
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt);

  CUDA_CALLABLE_MEMBER
  void swapP(int indx, int n);
  
  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP,std::size_t nStreams,libconfig::Config &cfg );
};

#endif
