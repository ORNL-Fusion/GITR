#ifndef _RECOMBINE_
#define _RECOMBINE_

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

struct recombine { 
	           const double dt;
        recombine(double _dt) : dt(_dt) {}
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(Particle &p) const { 
	if(p.hitWall == 0.0){        
	double trec;
	double Prec;
	double CoeffsRe[10] = {6.1164e-18, 5.1210e-17,1.2572e-16,1.9488e-16,2.5243e-16,3.0015e-16,3.3111e-16,3.5005e-16,3.5897e-16,3.5928e-16};
	double density = 1e19;
	double Temp_eV = 20;
	
	if(p.Z > 0)
	{
	trec = 1/(CoeffsRe[int(floor(p.Z- 0.5f))]*density);
	Prec = 1-exp(-dt/trec);
	}
	else 
	{
	Prec = 0.0;
	}
	
#if PARTICLESEEDS > 0	
	#ifdef __CUDACC__
	double r2 = curand_uniform(&p.streams[1]);
	#else
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double r2=dist(p.streams[1]);
	#endif
#else
    double r2 = 0.0;
#endif
						
	if(r2 <= Prec)
	{
		p.Z = p.Z-1;
	}         
	}	

	} 
};

#endif
