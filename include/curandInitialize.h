#ifndef _CURANDINITIAL_
#define _CURANDINITIAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif
#include <cstdlib>
#include <iostream>
#include "Particles.h"
#include "libconfig.h++"

template <typename T=std::mt19937>
struct curandInitialize{
#if USE_CUDA > 0
  curandState *s;
#else
  T *s;
#endif
  //int seed;
  bool fixed;
  
  curandInitialize(
#if USE_CUDA > 0
    curandState *_s,
#else
    T *_s,
#endif
    bool _fixed) : s(_s), fixed(_fixed) {}
    
  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx)
  {

#if USE_CUDA > 0

    if( fixed )
    {
      curand_init(indx, 0, 0, &s[indx]);
    }
    else
    {
      curand_init( clock64(), 0, 0, &s[indx] );
    }

#else

    if (fixed)
    {
      T s0(1234^indx);
      s[indx] = s0;
    }

    else
    {
      std::random_device randDevice;
      std::mt19937 s0(randDevice());
      s[indx] = s0;
    }

#endif
  }
};
#endif
