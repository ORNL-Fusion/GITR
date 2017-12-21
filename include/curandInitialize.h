#ifndef _CURANDINITIAL_
#define _CURANDINITIAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include "Particle.h"
#include "libconfig.h++"


#if __CUDACC__
//__global__
struct curandInitialize{
   curandState *s;
  int seed;
  
 curandInitialize(curandState *_s, int _seed) : s(_s), seed(_seed) {} 
    CUDA_CALLABLE_MEMBER_DEVICE
//void curandInitialize(curandState *s, int seed){
    void operator()(std::size_t indx) const {

       uint32_t block_id=blockIdx.y*gridDim.x+blockIdx.x;
       uint32_t blockSize=blockDim.z*blockDim.y*blockDim.x;
       uint32_t thread_id=threadIdx.z*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x;
       int seedTotal;
       int deviceN;
      // cudaGetDevice(&deviceN);
       //if (thread_id > 639){
       //    seedTotal = thread_id*2;
       //}
       //else{ seedTotal = thread_id;}

        curand_init(indx, 0, 0, &s[indx]);
    }
    };
#endif
#endif
