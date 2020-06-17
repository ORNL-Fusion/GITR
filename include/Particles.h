#ifndef _PARTICLES_
#define _PARTICLES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include "flags.hpp"
#include <libconfig.h++>

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

class Particles : public ManagedAllocation {
public:
  std::size_t nParticles;
  sim::Array<int> index;
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
  sim::Array<float> nu_s;
  sim::Array<float> vD;
  sim::Array<int> tt;
  sim::Array<int> hasLeaked;
  sim::Array<float> leakZ;
#ifdef __CUDACC__
  sim::Array<curandState> stream_ionize;
  //sim::Array<curandState> streams_rec;
  //sim::Array<curandState> streams_collision1;
  //sim::Array<curandState> streams_collision2;
  //sim::Array<curandState> streams_collision3;
  //sim::Array<curandState> streams_diff;
  //sim::Array<curandState> streams_surf;
#else
  sim::Array<std::mt19937> stream_ionize;
  //sim::Array<std::mt19937> streams_rec;
  //sim::Array<std::mt19937> streams_collision1;
  //sim::Array<std::mt19937> streams_collision2;
  //sim::Array<std::mt19937> streams_collision3;
  //sim::Array<std::mt19937> streams_diff;
  //sim::Array<std::mt19937> streams_surf;
#endif

  sim::Array<float> hitWall;
  sim::Array<int> surfaceHit;
  sim::Array<int> firstCollision;
  sim::Array<float> transitTime;
  sim::Array<float> distTraveled;
  sim::Array<int> wallIndex;
  sim::Array<float> perpDistanceToSurface;
  sim::Array<float> test;
  sim::Array<float> test0;
  sim::Array<float> test1;
  sim::Array<float> test2;
  sim::Array<float> test3;
  sim::Array<float> test4;
  sim::Array<float> distanceTraveled;
  sim::Array<float> weight;
  sim::Array<float> firstIonizationZ;
  sim::Array<float> firstIonizationT;

  //  void BorisMove(double dt, double xMinV, double xMaxV, double yMin, double
  //  yMax, double zMin, double zMax);

  //  void Ionization(double dt);
  CUDA_CALLABLE_MEMBER
  void setParticle(int indx, float x, float y, float z, float Ex, float Ey,
                   float Ez, float Z, float amu, float charge) {

    // this->index[indx] = indx;
    this->xprevious[indx] = x;
    this->yprevious[indx] = y;
    this->zprevious[indx] = z;
    this->x[indx] = x;
    this->y[indx] = y;
    this->z[indx] = z;
    this->Z[indx] = Z;
    this->charge[indx] = charge;
    this->amu[indx] = amu;
    this->hitWall[indx] = 0.0;
    this->wallIndex[indx] = 0;
    //        float Ex,Ey,Ez;
    //        Ex = E*cos(theta)*sin(phi);
    //        Ey = E*sin(theta)*sin(phi);
    //        Ez = E*cos(phi);
    this->vx[indx] =
        Ex / std::abs(Ex) *
        std::sqrt(2.0 * std::abs(Ex) * 1.60217662e-19 / (amu * 1.6737236e-27));
    this->vy[indx] =
        Ey / std::abs(Ey) *
        std::sqrt(2.0 * std::abs(Ey) * 1.60217662e-19 / (amu * 1.6737236e-27));
    this->vz[indx] =
        Ez / std::abs(Ez) *
        std::sqrt(2.0 * std::abs(Ez) * 1.60217662e-19 / (amu * 1.6737236e-27));

    if (Ex == 0.0)
      this->vx[indx] = 0.0;
    if (Ey == 0.0)
      this->vy[indx] = 0.0;
    if (Ez == 0.0)
      this->vz[indx] = 0.0;
    // std::cout << " velocity " << this->vx[indx] << " " << this->vy[indx] << "
    // " << this->vz[indx] << std::endl;
  };

  CUDA_CALLABLE_MEMBER
  void setParticleV(int indx, float x, float y, float z, float Vx, float Vy,
                    float Vz, float Z, float amu, float charge) {
    int indTmp = indx;
    this->index[indx] = indTmp;
    this->xprevious[indx] = x;
    this->yprevious[indx] = y;
    this->zprevious[indx] = z;
    this->x[indx] = x;
    this->y[indx] = y;
    this->z[indx] = z;
    this->Z[indx] = Z;
    this->charge[indx] = charge;
    this->amu[indx] = amu;
    this->hitWall[indx] = 0.0;
    this->wallIndex[indx] = 0;
    this->vx[indx] = Vx;
    this->vy[indx] = Vy;
    this->vz[indx] = Vz;
    this->v[indx] = std::sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
  };
  CUDA_CALLABLE_MEMBER
  void swapP(int indx, int n) {
    int iT = this->index[indx];
    float xpT = this->xprevious[indx];
    float ypT = this->yprevious[indx];
    float zpT = this->zprevious[indx];
    float xT = this->x[indx];
    float yT = this->y[indx];
    float zT = this->z[indx];
    float wT = this->weight[indx];
    float ZT = this->Z[indx];
    float cT = this->charge[indx];
    float aT = this->amu[indx];
    float hWT = this->hitWall[indx];
    int wIT = this->wallIndex[indx];
    float vxT = this->vx[indx];
    float vyT = this->vy[indx];
    float vzT = this->vz[indx];
    int wHT = this->surfaceHit[indx];
    float ttT = this->transitTime[indx];
    float dtT = this->distTraveled[indx];
    float firstIonizationZT = this->firstIonizationZ[indx];
    float firstIonizationTT = this->firstIonizationT[indx];

    this->index[indx] = this->index[n];
    this->xprevious[indx] = this->xprevious[n];
    this->yprevious[indx] = this->yprevious[n];
    this->zprevious[indx] = this->zprevious[n];
    this->x[indx] = this->x[n];
    this->y[indx] = this->y[n];
    this->z[indx] = this->z[n];
    this->weight[indx] = this->weight[n];
    this->Z[indx] = this->Z[n];
    this->charge[indx] = this->charge[n];
    this->amu[indx] = this->amu[n];
    this->hitWall[indx] = this->hitWall[n];
    this->wallIndex[indx] = this->wallIndex[n];
    this->vx[indx] = this->vx[n];
    this->vy[indx] = this->vy[n];
    this->vz[indx] = this->vz[n];
    this->surfaceHit[indx] = this->surfaceHit[n];
    this->transitTime[indx] = this->transitTime[n];
    this->distTraveled[indx] = this->distTraveled[n];
    this->firstIonizationZ[indx] = this->firstIonizationZ[n];
    this->firstIonizationT[indx] = this->firstIonizationT[n];

    this->index[n] = iT;
    this->xprevious[n] = xpT;
    this->yprevious[n] = ypT;
    this->zprevious[n] = zpT;
    this->x[n] = xT;
    this->y[n] = yT;
    this->z[n] = zT;
    this->weight[n] = wT;
    this->Z[n] = ZT;
    this->charge[n] = cT;
    this->amu[n] = aT;
    this->hitWall[n] = hWT;
    this->wallIndex[n] = wIT;
    this->vx[n] = vxT;
    this->vy[n] = vyT;
    this->vz[n] = vzT;
    this->surfaceHit[n] = wHT;
    this->transitTime[n] = ttT;
    this->distTraveled[n] = dtT;
    this->firstIonizationZ[n] = firstIonizationZT;
    this->firstIonizationT[n] = firstIonizationTT;
  };
  
  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP,std::size_t nStreams,libconfig::Config &cfg, Flags *gitr_flags) : 
    nParticles{getVariable_cfg<unsigned int> (cfg,"impurityParticleSource.nP")},
    index{nParticles, 0}, 
    x{nP,0.0}, 
    y{nP,0.0}, 
    z{nP,0.0}, 
    xprevious{nParticles,0.0},
    yprevious{nParticles,0.0}, 
    zprevious{nParticles,0.0},
    v{nParticles, 0.0}, 
    vx{nParticles,std::sqrt(static_cast<float>(2.0*getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.energy_eV")*
        1.602e-19/getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
    std::cos(getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.theta"))*
    std::sin(getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.phi"))}, 
    vy{nParticles,std::sqrt(static_cast<float>(2.0*getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.energy_eV")*
		    1.602e-19/getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
    std::sin(getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.theta"))*
    std::sin(getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.phi"))}, 
    vz{nParticles,std::sqrt(static_cast<float>(2.0*getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.energy_eV")*
		    1.602e-19/getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
    std::cos(getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.phi"))}, 
    Z{nParticles},
    amu{nParticles,getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.impurity_amu")}, 
    charge{nParticles,getVariable_cfg<float> (cfg,"impurityParticleSource.initialConditions.charge")}, 
    newVelocity{nParticles}, 
    nu_s{nParticles}, 
    vD{nParticles, 0.0}, 
    tt{nParticles, 0}, 
    hasLeaked{nParticles, 0}, 
    leakZ{nParticles,0.0}, 
    stream_ionize{nParticles},

//streams_rec{nStreams},streams_collision1{nStreams},streams_collision2{nStreams},
  //    streams_collision3{nStreams},streams_diff{nStreams},streams_surf{nStreams},
    hitWall{nParticles, 0.0},
    surfaceHit{nParticles, -1},
    firstCollision{nParticles, 1}, 
    transitTime{nParticles, 0.0}, 
    distTraveled{nParticles, 0.0},
    wallIndex{nParticles,0},
    perpDistanceToSurface{nParticles,0.0}, 
    test{nParticles, 0.0}, 
    test0{nParticles, 0.0},
    test1{nParticles, 0.0}, 
    test2{nParticles, 0.0}, 
    test3{nParticles, 0.0}, 
    test4{nParticles, 0.0},
    distanceTraveled{nParticles,0.0}, 
    weight{nParticles, 1.0}, 
    firstIonizationZ{nParticles, 0.0},
    firstIonizationT{nParticles, 0.0} {};
  
  //CUDA_CALLABLE_MEMBER_DEVICE
  //      #if __CUDACC__  
  //sim::Array<curandState> initialize_random_streams(int nStreams,libconfig::Config &cfg, Flags *gitr_flags)
  //{
  //        sim::Array<curandState> stream(nStreams);	  
  //      #else      
  //sim::Array<std::mt19937> initialize_random_streams(int nStreams,libconfig::Config &cfg, Flags *gitr_flags)
  //{
  //  sim::Array<std::mt19937> stream(nStreams);
  //#endif
//
//    if(gitr_flags->FIXED_SEEDS)
//    {
//      int seed0 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
//      for (int i =0;i < nParticles; i++) 
//      {
//        #if __CUDACC__  
//          curand_init(i, 0, 0, &stream[i]);
//        #else      
//          std::mt19937 s0(seed0+i);
//          stream_ionize[i] = s0;
//        #endif
//      }
//    }
//    else
//    {
//      
//      std::random_device randDeviceInit;
//      for (int i =0;i < nStreams; i++) 
//      {
//        #if __CUDACC__  
//          curand_init(clock64(), 0, 0, &stream[i]);
//        #else      
//          std::mt19937 s0(randDeviceInit());
//          stream[i] = s0;
//        #endif
//      }
//    }
//     
//    return stream;
//  }

};

#endif
