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
  sim::Array<float> vx;
  sim::Array<float> vy;
  sim::Array<float> vz;
  sim::Array<float> Z;
  sim::Array<float> amu;
  sim::Array<float> charge;
  sim::Array<float> newVelocity;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
  //sim::Array<curandState> streams;
  //sim::Array<curandState> streams_rec;
  //sim::Array<curandState> streams_collision1;
  //sim::Array<curandState> streams_collision2;
  //sim::Array<curandState> streams_collision3;
  //sim::Array<curandState> streams_diff;
  //sim::Array<curandState> streams_surf;
#else
  //sim::Array<std::mt19937> streams;
  //sim::Array<std::mt19937> streams_rec;
  //sim::Array<std::mt19937> streams_collision1;
  //sim::Array<std::mt19937> streams_collision2;
  //sim::Array<std::mt19937> streams_collision3;
  //sim::Array<std::mt19937> streams_diff;
  //sim::Array<std::mt19937> streams_surf;
#endif
#endif

  sim::Array<float> hitWall;
  sim::Array<int> wallHit;
  sim::Array<float> transitTime;
  sim::Array<int> wallIndex;
  sim::Array<float> perpDistanceToSurface;
  sim::Array<float> test;
  sim::Array<float> test0;
  sim::Array<float> test1;
  sim::Array<float> test2;
  sim::Array<float> test3;
  sim::Array<float> distanceTraveled;
  sim::Array<float> weight;
  sim::Array<float> PionizationPrevious;
  sim::Array<float> PrecombinationPrevious;
  sim::Array<float> firstIonizationZ;
  sim::Array<float> firstIonizationT;

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
//        float Ex,Ey,Ez;
//        Ex = E*cos(theta)*sin(phi);
//        Ey = E*sin(theta)*sin(phi);
//        Ez = E*cos(phi);
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

      };    
  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP) :
   nParticles{nP}, x{nP}, y{nP}, z{nP}, xprevious{nP}, yprevious{nP}, zprevious{nP},
   vx{nP}, vy{nP}, vz{nP}, Z{nP}, amu{nP}, charge{nP}, newVelocity{nP},
#if PARTICLESEEDS > 0
  //    streams{nP},streams_rec{nP},streams_collision1{nP},streams_collision2{nP},
  //    streams_collision3{nP},streams_diff{nP},streams_surf{nP},
#endif
      hitWall{nP,0.0},
   transitTime{nP,0.0},wallHit{nP,0}, wallIndex{nP}, perpDistanceToSurface{nP}, 
      test{nP},test0{nP},test1{nP},test2{nP},test3{nP},distanceTraveled{nP},weight{nP,1.0}, PionizationPrevious{nP,1.0},
    PrecombinationPrevious{nP,1.0}, firstIonizationZ{nP,0.0},firstIonizationT{nP,0.0} {};   

};

#endif
