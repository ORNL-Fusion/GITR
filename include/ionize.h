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
#include "Fields.h"
#include "flags.hpp"
#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/random.h>
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

#if __CUDA_ARCH__
CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(curandState *state,int indx)
{
        curandState localState = state[indx];
        float r1 = curand_uniform(&localState);
        state[indx] = localState;
	return r1;
}
#endif

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(std::mt19937 *state,int indx)
{
        std::uniform_real_distribution<float> dist(0.0, 1.0);
	std::mt19937 this_state = state[indx];
        float r1 = dist(this_state);
	state[indx] = this_state;
	return r1;
}

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(std::minstd_rand *state,int indx)
{
        std::uniform_real_distribution<float> dist(0.0, 1.0);
	std::minstd_rand this_state = state[indx];
        float r1 = dist(this_state);
	state[indx] = this_state;
	return r1;
}

#if __CUDA_ARCH__
CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(curandState *state,int indx)
{
        curandState localState = state[indx];
        double r1 = curand_uniform_double(&localState);
        state[indx] = localState;
	return r1;
}
#endif

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(std::mt19937 *state,int indx)
{
        std::uniform_real_distribution<double> dist(0.0, 1.0);
	std::mt19937 this_state = state[indx];
        double r1 = dist(this_state);
	state[indx] = this_state;
	return r1;
}

template <typename T=std::mt19937> 
struct ionize {
  Flags *flags;
  Particles *particlesPointer;
  int nR_Dens;
  int nZ_Dens;
  gitr_precision *DensGridr;
  gitr_precision *DensGridz;
  gitr_precision *ne;
  int nR_Temp;
  int nZ_Temp;
  gitr_precision *TempGridr;
  gitr_precision *TempGridz;
  gitr_precision *te;
  int nTemperaturesIonize;
  int nDensitiesIonize;
  gitr_precision *gridDensity_Ionization;
  gitr_precision *gridTemperature_Ionization;
  gitr_precision *rateCoeff_Ionization;
  gitr_precision dt;
  gitr_precision tion;
  void (ionize::*func)(std::size_t);
  // int& tt;
  int xx1;
  T *state;
  //Field *field1;

  gitr_precision  * random_uniform_number;
//#if __CUDA_ARCH__
//  curandState *state;
//#else
//  std::mt19937 *state;
//#endif
  //ionize() : dt(0.0){};
  ionize(Flags *_flags, Particles *_particlesPointer, gitr_precision _dt,T *_state,
//#if __CUDA_ARCH__
//         curandState *_state,
//#else
//         std::mt19937 *_state,
//#endif
         int _nR_Dens, int _nZ_Dens, gitr_precision *_DensGridr, gitr_precision *_DensGridz,
         gitr_precision *_ne, int _nR_Temp, int _nZ_Temp, gitr_precision *_TempGridr,
         gitr_precision *_TempGridz, gitr_precision *_te, int _nTemperaturesIonize,
         int _nDensitiesIonize, gitr_precision *_gridTemperature_Ionization,
         gitr_precision *_gridDensity_Ionization, gitr_precision *_rateCoeff_Ionization, gitr_precision *  _random_uniform_number)
      :

        flags(_flags), particlesPointer(_particlesPointer), nR_Dens(_nR_Dens),
        nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),
        ne(_ne), nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr),
        TempGridz(_TempGridz), te(_te),
        nTemperaturesIonize(_nTemperaturesIonize),
        nDensitiesIonize(_nDensitiesIonize),
        gridDensity_Ionization(_gridDensity_Ionization),
        gridTemperature_Ionization(_gridTemperature_Ionization),
        rateCoeff_Ionization(_rateCoeff_Ionization),
        dt(_dt),
        state(_state),random_uniform_number{_random_uniform_number} {
  }

  CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx) {
      //std::cout << "index " <<indx  << std::endl;
    if (flags->USE_IONIZATION) {
      //std::cout << " charge xyz nR nZ " << particlesPointer->charge[indx] << 
      // " " <<  particlesPointer->x[indx] << " " << 
      //    particlesPointer->y[indx] << " " <<  particlesPointer->z[indx] << 
      //   " " <<  nR_Temp << " " << 
      //    nZ_Temp << std::endl;
          // TempGridr, TempGridz, te, DensGridr, DensGridz, ne,
          //nTemperaturesIonize, nDensitiesIonize, gridTemperature_Ionization,
          //gridDensity_Ionization, rateCoeff_Ionization);
    if (flags->USE_ADAPTIVE_DT) {
	    dt = particlesPointer->dt[indx];
    }
      tion = interpRateCoeff2d(
          particlesPointer->charge[indx], particlesPointer->x[indx],
          particlesPointer->y[indx], particlesPointer->z[indx], nR_Temp,
          nZ_Temp, TempGridr, TempGridz, te, DensGridr, DensGridz, ne,
          nTemperaturesIonize, nDensitiesIonize, gridTemperature_Ionization,
          gridDensity_Ionization, rateCoeff_Ionization);
      //gitr_precision interp1 = field1->interpolate(1.0,2.0,3.0);
      gitr_precision P = expf(-dt / tion);
      gitr_precision P1 = 1.0 - P;
      gitr_precision r1 = get_rand_double(state,indx);
      //printf("dt P1 r1 %f %f %f \n", dt, P1, r1);
    if (flags->USE_ADAPTIVE_DT) {
      if (particlesPointer->hitWall[indx] == 0.0 && particlesPointer->advance[indx]) {
        gitr_precision r1 = get_rand_double(state,indx);
	//random_uniform_number[0] = r1;
        if (r1 <= P1) {
          particlesPointer->charge[indx] = particlesPointer->charge[indx] + 1;
        }
      }
    }
	    else{
      if (particlesPointer->hitWall[indx] == 0.0) {
        gitr_precision r1 = get_rand_double(state,indx);
	//random_uniform_number[0] = r1;
        if (r1 <= P1) {
      //printf("did ionize, adding one \n");
          particlesPointer->charge[indx] = particlesPointer->charge[indx] + 1;
        }
      }
	    }
    }
  }
};

//template <typename T=std::mt19937> 
//__global__ void ionize_kernel(ionize<curandState> i0, std::size_t indx) {
//   i0(indx);
//}
#endif
