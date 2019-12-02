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
    float tion;
    void (ionize::*func)(std::size_t);
    //int& tt;
    int xx1;
#if __CUDACC__
    curandState *state;
#else
    std::mt19937 *state;
#endif
        ionize():dt(0.0) {};
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
    float* _rateCoeff_Ionization
              ) : 
   
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
                                         nTemperaturesIonize(_nTemperaturesIonize),
                                         nDensitiesIonize(_nDensitiesIonize),
                                         gridDensity_Ionization(_gridDensity_Ionization),
                                         gridTemperature_Ionization(_gridTemperature_Ionization),
                                         rateCoeff_Ionization(_rateCoeff_Ionization),
                                         dt(_dt), // JDL missing tion here?
                                         state(_state) {
  }
        CUDA_CALLABLE_MEMBER_DEVICE
        void operator1(std::size_t indx) {}
        CUDA_CALLABLE_MEMBER_DEVICE
	  void doWork(std::size_t indx) {
    // ‘*this’ capture mode tells compiler to make a copy
    // of the object
    auto lam1 = [=, *this] __device__ { return xx1+1; };
    //auto lam1 = [=, *this] __device__ (std::size_t indx) { operator(indx); };
	  }
        CUDA_CALLABLE_MEMBER_DEVICE 
        virtual void operator()(std::size_t indx)  { 
	//if(particlesPointer->hitWall[indx] == 0.0){        
        //std::cout << "interpolating rate coeff at "<< particlesPointer->x[indx] << " " << particlesPointer->y[indx] << " " << particlesPointer->z[indx] << std::endl;
       tion = interpRateCoeff2d ( particlesPointer->charge[indx], particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx],nR_Temp,nZ_Temp, TempGridr,TempGridz,te,DensGridr,DensGridz, ne,nTemperaturesIonize,nDensitiesIonize,gridTemperature_Ionization,gridDensity_Ionization,rateCoeff_Ionization );	
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
       //particlesPointer->test0[indx] = P; 
       //particlesPointer->test1[indx] = P1; 
       //particlesPointer->test2[indx] = r1; 
	    if(r1 <= P1)
	    {
		  particlesPointer->charge[indx] = particlesPointer->charge[indx]+1;}
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
};
struct ionizeNothing : ionize
{
	int member;
	ionizeNothing():member(1){}
        CUDA_CALLABLE_MEMBER_DEVICE 
        void operator()(std::size_t indx)  {} 
};	
#endif
