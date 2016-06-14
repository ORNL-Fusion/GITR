#ifndef _IONIZE_
#define _IONIZE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particle.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

struct ionize { 
	           const double dt;
        ionize(double _dt) : dt(_dt) {}
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(Particle &p) const { 
	if(p.hitWall == 0.0){        
	double tion;
	double P1;
	double Coeffs[10] = {3.0875e-13, 9.6970e-14,3.8631e-14,2.0649e-14,3.5021e-15,1.6037e-15,7.0230e-17,1.7442e-17,6.1966e-18,1.8790e-18};
	double density = 1e19;
	double Temp_eV = 20;
	tion = 1/(Coeffs[int(floor(p.Z+ 0.5f))]*density);
	
	P1 = 1-exp(-dt/tion);
	
	#ifdef __CUDACC__
	double r1 = curand_uniform(&p.streams[0]);
	#else
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double r1=dist(p.streams[0]);
	#endif

	if(r1 <= P1)
	{
		p.Z = p.Z+1;
	} 
	}	

	} 
};

#endif
