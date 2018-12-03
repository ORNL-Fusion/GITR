#ifndef _OMPPRINT_
#define _OMPPRINT_

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

struct ompPrint { 

        ompPrint(void) {} 
    
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(std::size_t indx) const { 
        
#if USE_CUDA
#else  
      int tid = omp_get_thread_num();
      std::cout << " thread number " <<tid << std::endl; 
		
#endif
	} 
};

#endif
