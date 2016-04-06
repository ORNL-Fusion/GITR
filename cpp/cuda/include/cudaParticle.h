#ifndef _CUDAPARTICLE_
#define _CUDAPARTICLE_

#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <cstdlib>
#include <cmath>
#include <thrust/random.h>
#include <curand_kernel.h>

class cudaParticle {
	public:
		float x;
		float y;
		float z;
      	float vx;
      	float vy;
      	float vz;
      	float Z;
      	float amu;
	curandState s;
	//thrust::default_random_engine stream;

		__host__ __device__ 
        cudaParticle() {
            x=0.0;
	        y=0.0;
	        z=0.0;
        };
		__host__ __device__
        cudaParticle(float x,float y, float z, float Ex, float Ey, float Ez, float Z, float amu){
    
            this->x = x;
            this->y = y;
            this->z = z;

	        this->Z = Z;
		    this->amu = amu;

		    this->vx = Ex/fabs(Ex)*sqrt(2.0*fabs(Ex)*1.60217662e-19/(amu*1.6737236e-27));
		    this->vy = Ey/fabs(Ey)*sqrt(2.0*fabs(Ey)*1.60217662e-19/(amu*1.6737236e-27));
		    this->vz = Ez/fabs(Ez)*sqrt(2.0*fabs(Ez)*1.60217662e-19/(amu*1.6737236e-27));
		    
		    if(Ex == 0.0) this->vx = 0.0;
		    if(Ey == 0.0) this->vy = 0.0;
		    if(Ez == 0.0) this->vz = 0.0;
	
	//this->stream = eng; 
        };
};

#endif

