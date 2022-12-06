#ifndef _RECOMBINE_
#define _RECOMBINE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include "ionize.h"
#include "flags.hpp"
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

#include "interpRateCoeff.hpp"
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

template <typename T=std::mt19937>
struct recombine {
  Particles *particlesPointer;
  int nR_Dens;
  int nZ_Dens;
  gitr_precision* DensGridr;
  gitr_precision* DensGridz;
  gitr_precision* ne;
  int nR_Temp;
  int nZ_Temp;
  gitr_precision* TempGridr;
  gitr_precision* TempGridz;
  gitr_precision* te;
  int nTemperaturesRecomb;
  int nDensitiesRecomb;
  gitr_precision* gridDensity_Recombination;
  gitr_precision* gridTemperature_Recombination;
  gitr_precision* rateCoeff_Recombination;
  Flags* gitr_flags;
  gitr_precision dt;
  gitr_precision tion;
  T *state;
  int cylsymm;

  recombine(Particles *_particlesPointer, gitr_precision _dt,
      T *_state,
     int _nR_Dens,int _nZ_Dens,gitr_precision* _DensGridr,
     gitr_precision* _DensGridz,gitr_precision* _ne,int _nR_Temp, int _nZ_Temp,
     gitr_precision* _TempGridr, gitr_precision* _TempGridz,gitr_precision* _te,int _nTemperaturesRecomb,
     int _nDensitiesRecomb,gitr_precision* _gridTemperature_Recombination,gitr_precision* _gridDensity_Recombination,
     gitr_precision* _rateCoeff_Recombination, Flags* _gitr_flags, int cylsymm_ ) : 
     particlesPointer(_particlesPointer),

     nR_Dens(_nR_Dens),
     nZ_Dens(_nZ_Dens),
     DensGridr(_DensGridr),
     DensGridz(_DensGridz),
     ne(_ne),
     nR_Temp(_nR_Temp),
     nZ_Temp(_nZ_Temp),
     TempGridr(_TempGridr),
     TempGridz(_TempGridz),
     te(_te),
     nTemperaturesRecomb(_nTemperaturesRecomb),
     nDensitiesRecomb(_nDensitiesRecomb),
     gridDensity_Recombination(_gridDensity_Recombination),
     gridTemperature_Recombination(_gridTemperature_Recombination),
     rateCoeff_Recombination(_rateCoeff_Recombination),
     gitr_flags(_gitr_flags),
     dt(_dt),
     state(_state),
     cylsymm( cylsymm_ ) {}
 
  
  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx)
  {
    gitr_precision P1 = 0.0;
    gitr_precision r1 = 1.0;
      
    if (gitr_flags->USE_ADAPTIVE_DT)
    {
      dt = particlesPointer->dt[indx];
    }
    
    if(particlesPointer->charge[indx] > 0)
    {
      tion = interpRateCoeff2d ( particlesPointer->charge[indx]-1, particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],nR_Temp,nZ_Temp, TempGridr,TempGridz,te,DensGridr,DensGridz, ne,nTemperaturesRecomb,nDensitiesRecomb,gridTemperature_Recombination,gridDensity_Recombination,rateCoeff_Recombination, cylsymm );
      gitr_precision P = exp(-dt/tion);
      P1 = 1.0-P;
      r1 = get_rand_double(state,indx);
    }

    if (gitr_flags->USE_ADAPTIVE_DT)
    {
      if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->advance[indx])
      {
        if(r1 <= P1)
        {
          particlesPointer->charge[indx] = particlesPointer->charge[indx]-1;
        }
      }
    }
    else
    {
      if(particlesPointer->hitWall[indx] == 0.0)
      {
        if(r1 <= P1)
        {
          particlesPointer->charge[indx] = particlesPointer->charge[indx]-1;
        }
      }
    }

  }
};

#endif
