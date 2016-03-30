#ifndef _CUDAPARTICLE_
#define _CUDAPARTICLE_

#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <cstdlib>

class cudaParticle {
	public:
		double x;
		float y;
		float z;
      		float vx;
      		float vy;
      		float vz;
      		float Z;
      		float amu;

		__host__ __device__ 
        cudaParticle() {
            x=0.0;
	    y=0.0;
	    z=0.0;
        };
};

#endif
