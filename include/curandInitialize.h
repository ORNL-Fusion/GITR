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
__global__
void curandInitialize(curandState *s, int seed)
    {
        curand_init(seed, 0, 0, s);
    };
#endif
#endif
