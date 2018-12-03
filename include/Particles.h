#ifndef _PARTICLES_
#define _PARTICLES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>
#include "array.h"
#include "managed_allocation.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#include <random>

//CUDA_CALLABLE_MEMBER

class Particles : public ManagedAllocation {
public: 
  std::size_t nParticles;  
  sim::Array<float> x;
  sim::Array<float> y;
  sim::Array<float> z;
  sim::Array<float> xprevious;
  sim::Array<float> yprevious;
  sim::Array<float> zprevious;
  sim::Array<float> v;
  sim::Array<float> vx;
  sim::Array<float> vy;
  sim::Array<float> vz;
  sim::Array<float> Z;
  sim::Array<float> amu;
  sim::Array<float> charge;
  sim::Array<float> newVelocity;
  sim::Array<float> nu_E0;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
  sim::Array<curandState> streams;
#else
  sim::Array<std::mt19937> streams;
#endif
#endif

  sim::Array<float> hitWall;
  sim::Array<int> wallHit;
  sim::Array<int> firstCollision;
  sim::Array<float> transitTime;
  sim::Array<int> wallIndex;
  sim::Array<float> perpDistanceToSurface;
  sim::Array<float> test;

//  void BorisMove(double dt, double xMinV, double xMaxV, double yMin, double yMax, double zMin, double zMax);

//  void Ionization(double dt);
  CUDA_CALLABLE_MEMBER
  void  setParticle(int indx, float x, float y, float z, float Ex, float Ey, float Ez, float Z, float amu, float charge) {

        this->xprevious[indx] = x;
        this->yprevious[indx] = y;
        this->zprevious[indx] = z;
        this->x[indx] = x;
        this->y[indx] = y;
        this->z[indx] = z;
        this->Z[indx] = Z;
        this->charge[indx]= charge;
        this->amu[indx] = amu;
        this->hitWall[indx] = 0.0;
        this->wallIndex[indx] = 0;
        this->vx[indx] = Ex / fabs(Ex) * sqrt(2.0 * fabs(Ex) * 1.60217662e-19 / (amu * 1.6737236e-27));
        this->vy[indx] = Ey / fabs(Ey) * sqrt(2.0 * fabs(Ey) * 1.60217662e-19 / (amu * 1.6737236e-27));
        this->vz[indx] = Ez / fabs(Ez) * sqrt(2.0 * fabs(Ez) * 1.60217662e-19 / (amu * 1.6737236e-27));

        if (Ex == 0.0) this->vx[indx] = 0.0;
        if (Ey == 0.0) this->vy[indx] = 0.0;
        if (Ez == 0.0) this->vz[indx] = 0.0;
        //std::cout << " velocity " << this->vx[indx] << " " << this->vy[indx] << " " << this->vz[indx] << std::endl;
     };    

  CUDA_CALLABLE_MEMBER
  void  setParticleV(int indx, float x, float y, float z, float Vx, float Vy, float Vz, float Z, float amu, float charge) {
        int indTmp=indx;
        this->index[indx] = indTmp;
        this->xprevious[indx] = x;
        this->yprevious[indx] = y;
        this->zprevious[indx] = z;
        this->x[indx] = x;
        this->y[indx] = y;
        this->z[indx] = z;
        this->Z[indx] = Z;
        this->charge[indx]= charge;
        this->amu[indx] = amu;
        this->hitWall[indx] = 0.0;
        this->wallIndex[indx] = 0;
        this->vx[indx] = Vx;
        this->vy[indx] = Vy;
        this->vz[indx] = Vz;
        this->v[indx] = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

      };    

  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP) :
   nParticles{nP}, index{nP,0}, x{nP}, y{nP}, z{nP}, xprevious{nP}, yprevious{nP}, zprevious{nP},
   vx{nP}, vy{nP}, vz{nP},v{nP,0.0},Z{nP}, amu{nP}, charge{nP}, newVelocity{nP},nu_E0{nP},
#if PARTICLESEEDS > 0
      streams{nP},
#endif
      hitWall{nP,0.0},
   transitTime{nP,0.0},distTraveled{nP,0.0},
      wallHit{nP,0},firstCollision{nP,1}, wallIndex{nP}, perpDistanceToSurface{nP}, 
      test{nP,0.0},test0{nP,0.0},test1{nP,0.0},test2{nP,0.0},test3{nP,0.0},test4{nP,0.0},distanceTraveled{nP},weight{nP,1.0}, PionizationPrevious{nP,1.0},
    PrecombinationPrevious{nP,1.0}, firstIonizationZ{nP,0.0},firstIonizationT{nP,0.0} {};   

};

#endif