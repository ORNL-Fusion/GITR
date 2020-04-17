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
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

#include "interpRateCoeff.hpp"

template <typename T=std::mt19937> 
struct recombine { 
  Particles *particlesPointer;
  int nR_Dens;
  int nZ_Dens;
  float* DensGridr;
  float* DensGridz;
  float* ne;
  int nR_Temp;
  int nZ_Temp;
  float* TempGridr;
  float* TempGridz;
  float* te;
  int nTemperaturesRecomb;
  int nDensitiesRecomb;
  float* gridDensity_Recombination;
  float* gridTemperature_Recombination;
  float* rateCoeff_Recombination;
  const float dt;
  float tion;
  T *state;

  recombine(Particles *_particlesPointer, float _dt,
      T *_state,
     int _nR_Dens,int _nZ_Dens,float* _DensGridr,
     float* _DensGridz,float* _ne,int _nR_Temp, int _nZ_Temp,
     float* _TempGridr, float* _TempGridz,float* _te,int _nTemperaturesRecomb,
     int _nDensitiesRecomb,float* _gridTemperature_Recombination,float* _gridDensity_Recombination,
     float* _rateCoeff_Recombination) : 
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
                                               dt(_dt), // JDL missing tion?
                                               state(_state) {
  }
 
  
  CUDA_CALLABLE_MEMBER_DEVICE 
  void operator()(std::size_t indx) { 
      //std::cout << "rec index " <<indx  << std::endl;
      //if(indx > 1000)
      //{
      // std::cout << "wtf " <<indx  << std::endl;
      //} 
  float P1 = 0.0f;
      if(particlesPointer->charge[indx] > 0)
    {
       tion = interpRateCoeff2d ( particlesPointer->charge[indx]-1, particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],nR_Temp,nZ_Temp, TempGridr,TempGridz,te,DensGridr,DensGridz, ne,nTemperaturesRecomb,nDensitiesRecomb,gridTemperature_Recombination,gridDensity_Recombination,rateCoeff_Recombination);
       float P = expf(-dt/tion);
       P1 = 1.0-P;
    }

  if(particlesPointer->hitWall[indx] == 0.0)
	{        
        float r1 = get_rand(state,indx);
	if(r1 <= P1)
	{
        particlesPointer->charge[indx] = particlesPointer->charge[indx]-1;
	}         
   }	

  } 
};

#endif
