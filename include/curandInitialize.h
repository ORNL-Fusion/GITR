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
#if __CUDACC__
   curandState *s;
#else
   T *s;
#endif
  int seed;
  
 curandInitialize(
#if __CUDACC__		 
         curandState *_s,
#else
         T *_s,
#endif
         int _seed) : s(_s), seed(_seed) {} 
//void curandInitialize(curandState *_s, int _seed) : s(_s),seed(_seed) {}
    CUDA_CALLABLE_MEMBER_DEVICE
    void operator()(std::size_t indx) {
#if __CUDACC__
       //uint32_t block_id=blockIdx.y*gridDim.x+blockIdx.x;
       //uint32_t blockSize=blockDim.z*blockDim.y*blockDim.x;
       //uint32_t thread_id=threadIdx.z*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x;
       //int seedTotal;
       //int deviceN;
      //// cudaGetDevice(&deviceN);
       ////if (thread_id > 639){
       ////    seedTotal = thread_id*2;
       ////}
       ////else{ seedTotal = thread_id;}

        curand_init(indx, 0, 0, &s[indx]);
#else
        std::random_device randDevice;
        //std::mt19937 s0(randDevice());
        T s0(1234^indx);
        s[indx] = s0;
#endif
    }
    };
#endif
