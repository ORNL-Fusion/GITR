#pragma once
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
float get_rand(curandState *state,int indx);

#else

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(std::mt19937 *state,int indx);

#if __CUDA_ARCH__
CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(curandState *state,int indx);
#endif

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(std::mt19937 *state,int indx);

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
         gitr_precision *_gridDensity_Ionization, gitr_precision *_rateCoeff_Ionization, gitr_precision *  _random_uniform_number);

  CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx);
};

//template <typename T=std::mt19937> 
//__global__ void ionize_kernel(ionize<curandState> i0, std::size_t indx) {
//   i0(indx);
//}
#endif
