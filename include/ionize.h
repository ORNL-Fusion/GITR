#ifndef _IONIZE_
#define _IONIZE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

#include "interpRateCoeff.hpp"

struct ionize { 
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
    int nTemperaturesIonize;
    int nDensitiesIonize;
    float* gridDensity_Ionization;
    float* gridTemperature_Ionization;
    float* rateCoeff_Ionization;
    const float dt;
    int tt;
#if __CUDACC__
    curandState *state;
#else
    std::mt19937 *state;
#endif

        ionize(Particles *_particlesPointer, float _dt,
#if __CUDACC__
                curandState *_state,
#else
                std::mt19937 *_state,
#endif
                int _nR_Dens,int _nZ_Dens,float* _DensGridr,
    float* _DensGridz,float* _ne,int _nR_Temp, int _nZ_Temp,
    float* _TempGridr, float* _TempGridz,float* _te,int _nTemperaturesIonize,
    int _nDensitiesIonize,float* _gridTemperature_Ionization,float* _gridDensity_Ionization,
    float* _rateCoeff_Ionization, int _tt
              ) : particlesPointer(_particlesPointer), dt(_dt), state(_state), nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),ne(_ne),
        nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz), te(_te),
   nTemperaturesIonize(_nTemperaturesIonize), nDensitiesIonize(_nDensitiesIonize), gridTemperature_Ionization(_gridTemperature_Ionization), gridDensity_Ionization(_gridDensity_Ionization), rateCoeff_Ionization(_rateCoeff_Ionization), tt(_tt) {}
    
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(std::size_t indx) const { 
	//if(particlesPointer->hitWall[indx] == 0.0){        
        //std::cout << "interpolating rate coeff at "<< particlesPointer->x[indx] << " " << particlesPointer->y[indx] << " " << particlesPointer->z[indx] << std::endl;
        float tion = interpRateCoeff2d ( particlesPointer->charge[indx], particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],nR_Temp,nZ_Temp, TempGridr,TempGridz,te,DensGridr,DensGridz, ne,nTemperaturesIonize,nDensitiesIonize,gridTemperature_Ionization,gridDensity_Ionization,rateCoeff_Ionization );	
    //float PiP = particlesPointer->PionizationPrevious[indx];
    float P = expf(-dt/tion);
    //particlesPointer->PionizationPrevious[indx] = PiP*P;
    float P1 = 1.0-P;
    //std::cout << "tion P P1 " << tion << " " << P << " " << P1 << " " << PiP<< std::endl;
    if(particlesPointer->hitWall[indx] == 0.0)
    {
        //std::cout << "calculating r1 " << std::endl;i
#if PARTICLESEEDS > 0
	#ifdef __CUDACC__
	  //float r1 = 0.5;//curand_uniform(&particlesPointer->streams[indx]);
      float r1 = curand_uniform(&state[indx]);
	#else
	std::uniform_real_distribution<float> dist(0.0, 1.0);
	float r1=dist(state[indx]);
	//particlesPointer->test[indx] = r1;
        //std::cout << " r1 " << r1 << std::endl;
    #endif
#else
  #if __CUDACC__
    curandState localState = state[thread_id];
    float r1 = curand_uniform(&localState);
    state[thread_id] = localState;
  #else
    std::uniform_real_distribution<float> dist(0.0, 1.0);
    float r1=dist(state[0]);
  #endif
#endif
    //if(tt == 722)
    //{
      //std::cout << "r1 " << r1 << " " << P1 << std::endl;
		//particlesPointer->charge[indx] = particlesPointer->charge[indx]+1;
       //particlesPointer->test[indx] = tion; 
       //particlesPointer->test0[indx] = P1; 
       //particlesPointer->test1[indx] = r1; 
	    if(r1 <= P1)
	    {
		  particlesPointer->charge[indx] = particlesPointer->charge[indx]+1;
          particlesPointer->PionizationPrevious[indx] = 1.0;
        //std::cout << "Particle " << indx << " ionized at step " << tt << std::endl;
       if(particlesPointer->firstIonizationZ[indx] == 0.0)
       {
           particlesPointer->firstIonizationZ[indx] = particlesPointer->z[indx];
       }
	    }
        else
        {
       if(particlesPointer->firstIonizationZ[indx] == 0.0)
       {
           particlesPointer->firstIonizationT[indx] = particlesPointer->firstIonizationT[indx] + dt;
       }

        } 
       
    //} 
	}	

	} 
};

#endif
