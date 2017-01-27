#ifndef _BOUNDARYMOD_
#define _BOUNDARYMOD_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>
#include "array.h"
#include "managed_allocation.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#include <random>

//CUDA_CALLABLE_MEMBER

class BoundaryModifiable : public ManagedAllocation {
public: 
  std::size_t nLines;  
  sim::Array<float> impacts;
  

  CUDA_CALLABLE_MEMBER
  BoundaryModifiable(std::size_t nLines) :
   impacts{nLines,0.0} {};   

};

#endif
