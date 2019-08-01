#ifndef _PARTICLES_
#define _PARTICLES_

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
  sim::Array<float> PionizationPrevious;
  sim::Array<float> PrecombinationPrevious;
  sim::Array<float> firstIonizationZ;
  sim::Array<float> firstIonizationT;

//  void BorisMove(double dt, double xMinV, double xMaxV, double yMin, double yMax, double zMin, double zMax);

//  void Ionization(double dt);
  CUDA_CALLABLE_MEMBER
  void  setParticle(int indx, float x, float y, float z, float Ex, float Ey, float Ez, float Z, float amu, float charge) {

        //this->index[indx] = indx;
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
  void  swapP(int indx,int n) {
        int iT=this->index[indx];
        float xpT=this->xprevious[indx];
        float ypT=this->yprevious[indx];
        float zpT=this->zprevious[indx];
        float xT=this->x[indx];
        float yT=this->y[indx];
        float zT=this->z[indx];
        float wT=this->weight[indx];
        float ZT=this->Z[indx];
        float cT=this->charge[indx];
        float aT=this->amu[indx];
        float hWT=this->hitWall[indx];
        int wIT=this->wallIndex[indx];
        float vxT=this->vx[indx];
        float vyT=this->vy[indx];
        float vzT=this->vz[indx];
        int wHT=this->wallHit[indx];
  float ttT=this->transitTime[indx];
  float dtT=this->distTraveled[indx];
  float firstIonizationZT=this->firstIonizationZ[indx];
  float firstIonizationTT=this->firstIonizationT[indx];

        this->index[indx] = this->index[n];
        this->xprevious[indx] = this->xprevious[n];
        this->yprevious[indx] = this->yprevious[n];
        this->zprevious[indx] = this->zprevious[n];
        this->x[indx] =this->x[n];
        this->y[indx] =this->y[n];
        this->z[indx] = this->z[n];
        this->weight[indx] = this->weight[n];
        this->Z[indx] = this->Z[n];
        this->charge[indx]= this->charge[n];
        this->amu[indx] =this->amu[n];
        this->hitWall[indx] =this->hitWall[n];
        this->wallIndex[indx]=this->wallIndex[n];
        this->vx[indx] =this->vx[n];
        this->vy[indx] = this->vy[n];
        this->vz[indx] =this->vz[n];
        this->wallHit[indx]=this->wallHit[n];
        this->transitTime[indx]=this->transitTime[n];
        this->distTraveled[indx]=this->distTraveled[n];
        this->firstIonizationZ[indx]=this->firstIonizationZ[n];
        this->firstIonizationT[indx]=this->firstIonizationT[n];

        this->index[n] = iT;
        this->xprevious[n] = xpT;
        this->yprevious[n] = ypT;
        this->zprevious[n] = zpT;
        this->x[n] =xT;
        this->y[n] =yT;
        this->z[n] =zT;
        this->weight[n] = wT;
        this->Z[n] = ZT;
        this->charge[n]= cT;
        this->amu[n] =aT;
        this->hitWall[n] =hWT;
        this->wallIndex[n]=wIT;
        this->vx[n] =vxT;
        this->vy[n] =vyT;
        this->vz[n] =vzT;
        this-> wallHit[n]=wHT;
        this->transitTime[n]=ttT;
        this->distTraveled[n]=dtT;
        this->firstIonizationZ[n]=firstIonizationZT;
        this->firstIonizationT[n]=firstIonizationTT;
      };    
  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP) :
   nParticles{nP}, index{nP,0}, x{nP}, y{nP}, z{nP}, xprevious{nP}, yprevious{nP}, zprevious{nP},
   vx{nP}, vy{nP}, vz{nP},v{nP,0.0},Z{nP}, amu{nP}, charge{nP}, newVelocity{nP},nu_s{nP},vD{nP,0.0},
#if PARTICLESEEDS > 0
  //    streams{nP},streams_rec{nP},streams_collision1{nP},streams_collision2{nP},
  //    streams_collision3{nP},streams_diff{nP},streams_surf{nP},
#endif
      hitWall{nP,0.0},
   transitTime{nP,0.0},distTraveled{nP,0.0},
      wallHit{nP,0},firstCollision{nP,1}, wallIndex{nP}, perpDistanceToSurface{nP}, 
      test{nP,0.0},test0{nP,0.0},test1{nP,0.0},test2{nP,0.0},test3{nP,0.0},test4{nP,0.0},distanceTraveled{nP},weight{nP,1.0}, PionizationPrevious{nP,1.0},
    PrecombinationPrevious{nP,1.0}, firstIonizationZ{nP,0.0},firstIonizationT{nP,0.0}, tt{nP,0},hasLeaked{nP,0},leakZ{nP,0.0} {};   

};

#endif
