#ifndef _PARTICLES_
#define _PARTICLES_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include "utils.h"
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
  sim::Array<gitr_precision> transitAngle;
  sim::Array<gitr_precision> distTraveled;
  sim::Array<int> wallIndex;
  sim::Array<gitr_precision> perpDistanceToSurface;
  sim::Array<gitr_precision> test;
  sim::Array<gitr_precision> test0;
  sim::Array<gitr_precision> test1;
  sim::Array<gitr_precision> test2;
  sim::Array<gitr_precision> angle;
  sim::Array<gitr_precision> max_z;
  sim::Array<gitr_precision> distanceTraveled;
  sim::Array<gitr_precision> weight;
  sim::Array<gitr_precision> firstIonizationZ;
  sim::Array<gitr_precision> firstIonizationT;
  sim::Array<gitr_precision> dt;
  sim::Array<gitr_precision> time;
  sim::Array<bool> advance;
  sim::Array<gitr_precision> f_psi;

  CUDA_CALLABLE_MEMBER
  void setParticle(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                   gitr_precision Ex, gitr_precision Ey, gitr_precision Ez,
                   gitr_precision Z, gitr_precision amu, gitr_precision charge) 
  {

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
  };

  CUDA_CALLABLE_MEMBER
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt)
  {
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
    this->dt[indx] = dt;
  };

  CUDA_CALLABLE_MEMBER
  void swapP(int indx, int n)
  {
    int iT = this->index[indx];
    gitr_precision xT = this->x[indx];
    gitr_precision yT = this->y[indx];
    gitr_precision zT = this->z[indx];
    gitr_precision xpT = this->xprevious[indx];
    gitr_precision ypT = this->yprevious[indx];
    gitr_precision zpT = this->zprevious[indx];
    gitr_precision vT = this->v[indx];
    gitr_precision vxT = this->vx[indx];
    gitr_precision vyT = this->vy[indx];
    gitr_precision vzT = this->vz[indx];
    gitr_precision ZT = this->Z[indx];
    gitr_precision cT = this->charge[indx];
    gitr_precision aT = this->amu[indx];
    gitr_precision nvT = this->newVelocity[indx];
    gitr_precision nsT = this->nu_s[indx];
    gitr_precision vdT = this->vD[indx];
    int hlT = this->hasLeaked[indx];
    gitr_precision lzT = this->leakZ[indx];
    gitr_precision max_zT = this->max_z[indx];
    gitr_precision angleT = this->angle[indx];
    gitr_precision wT = this->weight[indx];
    gitr_precision hWT = this->hitWall[indx];
    int wIT = this->wallIndex[indx];
    int wHT = this->surfaceHit[indx];
    int fcT = this->firstCollision[indx];
    gitr_precision ttT = this->transitTime[indx];
    gitr_precision ttA = this->transitAngle[indx];
    gitr_precision dtT = this->distTraveled[indx];
    gitr_precision firstIonizationZT = this->firstIonizationZ[indx];
    gitr_precision firstIonizationTT = this->firstIonizationT[indx];
    gitr_precision pdtsT = this->perpDistanceToSurface[indx];
    gitr_precision dtvT = this->distanceTraveled[indx];
    gitr_precision dt_T = this->dt[indx];
    gitr_precision timeT = this->time[indx];
    bool advanceT = this->advance[indx];

    this->index[indx] = this->index[n];
    this->x[indx] = this->x[n];
    this->y[indx] = this->y[n];
    this->z[indx] = this->z[n];
    this->xprevious[indx] = this->xprevious[n];
    this->yprevious[indx] = this->yprevious[n];
    this->zprevious[indx] = this->zprevious[n];
    this->v[indx] = this->v[n];
    this->vx[indx] = this->vx[n];
    this->vy[indx] = this->vy[n];
    this->vz[indx] = this->vz[n];
    this->Z[indx] = this->Z[n];
    this->charge[indx] = this->charge[n];
    this->amu[indx] = this->amu[n];
    this->newVelocity[indx] = this->newVelocity[n];
    this->nu_s[indx] = this->nu_s[n];
    this->vD[indx] = this->vD[n];
    this->hasLeaked[indx] = this->hasLeaked[n];
    this->leakZ[indx] = this->leakZ[n];
    this->max_z[indx] = this->max_z[n];
    this->angle[indx] = this->angle[n];
    this->weight[indx] = this->weight[n];
    this->hitWall[indx] = this->hitWall[n];
    this->wallIndex[indx] = this->wallIndex[n];
    this->surfaceHit[indx] = this->surfaceHit[n];
    this->firstCollision[indx] = this->firstCollision[n];
    this->transitTime[indx] = this->transitTime[n];
    this->transitAngle[indx] = this->transitAngle[n];
    this->distTraveled[indx] = this->distTraveled[n];
    this->firstIonizationZ[indx] = this->firstIonizationZ[n];
    this->firstIonizationT[indx] = this->firstIonizationT[n];
    this->dt[indx] = this->dt[n];
    this->time[indx] = this->time[n];
    this->advance[indx] = this->advance[n];

    this->index[n] = iT;
    this->x[n] = xT;
    this->y[n] = yT;
    this->z[n] = zT;
    this->xprevious[n] = xpT;
    this->yprevious[n] = ypT;
    this->zprevious[n] = zpT;
    this->v[n] = vT;
    this->vx[n] = vxT;
    this->vy[n] = vyT;
    this->vz[n] = vzT;
    this->Z[n] = ZT;
    this->charge[n] = cT;
    this->amu[n] = aT;
    this->newVelocity[n] = nvT;
    this->nu_s[n] = nsT;
    this->vD[n] = vdT;
    this->hasLeaked[n] = hlT;
    this->leakZ[n] = lzT;
    this->max_z[n] = max_zT;
    this->angle[n] = angleT;
    this->weight[n] = wT;
    this->hitWall[n] = hWT;
    this->wallIndex[n] = wIT;
    this->surfaceHit[n] = wHT;
    this->firstCollision[n] = fcT;
    this->transitTime[n] = ttT;
    this->transitAngle[n] = ttA;
    this->distTraveled[n] = dtT;
    this->firstIonizationZ[n] = firstIonizationZT;
    this->firstIonizationT[n] = firstIonizationTT;
    this->dt[n] = dt_T;
    this->time[n] = timeT;
    this->advance[n] = advanceT;
  };
  
  CUDA_CALLABLE_MEMBER
  Particles(std::size_t nP,std::size_t nStreams,libconfig::Config &cfg) :
    nParticles{getVariable_cfg<unsigned int> (cfg,"impurityParticleSource.nP")},
    index{nParticles, 0},
    x{nP,0.0},
    y{nP,0.0},
    z{nP,0.0},
    xprevious{nParticles,0.0},
    yprevious{nParticles,0.0},
    zprevious{nParticles,0.0},
    v{nParticles, 0.0},
    vx{nParticles,std::sqrt(static_cast<gitr_precision>(2.0*getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.energy_eV")*
          1.602e-19/getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
          std::cos(getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.theta"))*
          std::sin(getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.phi"))},
    vy{nParticles,std::sqrt(static_cast<gitr_precision>(2.0*getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.energy_eV")*
          1.602e-19/getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
          std::sin(getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.theta"))*
          std::sin(getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.phi"))},
    vz{nParticles,std::sqrt(static_cast<gitr_precision>(2.0*getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.energy_eV")*
          1.602e-19/getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.impurity_amu")/1.66e-27))*
          std::cos(getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.phi"))},
    Z{nParticles,0.0},
    amu{nParticles,getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.impurity_amu")},
    charge{nParticles,getVariable_cfg<gitr_precision>
          (cfg,"impurityParticleSource.initialConditions.charge")},
    newVelocity{nParticles},
    nu_s{nParticles},
    vD{nParticles, 0.0},
    tt{nParticles, 0},
    hasLeaked{nParticles, 0},
    leakZ{nParticles,0.0},
    stream_ionize{nParticles},
    hitWall{nParticles, 0.0},
    surfaceHit{nParticles, -1},
    firstCollision{nParticles, 1},
    transitTime{nParticles, 0.0},
    transitAngle{nParticles, 0.0},
    distTraveled{nParticles, 0.0},
    wallIndex{nParticles,0},
    perpDistanceToSurface{nParticles,0.0},
    test{nParticles, 0.0},
    test0{nParticles, 0.0},
    test1{nParticles, 0.0},
    test2{nParticles, 0.0},
    max_z{nParticles, 0.0},
    angle{nParticles, 0.0},
    distanceTraveled{nParticles,0.0},
    weight{nParticles, 1.0},
    firstIonizationZ{nParticles, 0.0},
    firstIonizationT{nParticles, 0.0},
    dt{nParticles,0.0},
    time{nParticles,0.0},
    advance{nParticles,false},
    f_psi{nParticles,1.0} {};
};

#endif
