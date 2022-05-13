#ifndef _PARTICLE_
#define _PARTICLE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <cmath>
#include <stdio.h>
//#include <vector>
#include "array.h"
//#include "managed_allocation.h"

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#else
#include <random>
#endif

//CUDA_CALLABLE_MEMBER

class Particle  {
	public:
	    float x;
	    float y;
	    float z;
        float xprevious;
        float yprevious;
        float zprevious;
      	float vx;
      	float vy;
      	float vz;
      	float Z;
      	float amu;
        float charge;
    #ifdef __CUDACC__
	curandState streams[7];
	#else
        std::mt19937 streams[7];
        #endif
	float hitWall;
    float transitTime;
    int wallIndex;
    float perpDistanceToSurface;
	float seed0;
	void BorisMove(double dt, double xMinV,double xMaxV,double yMin,double yMax,double zMin,double zMax);
	void Ionization(double dt);

	CUDA_CALLABLE_MEMBER
        Particle() {
            x=0.0;
	    y=0.0;
	    z=0.0;
        };
	
	CUDA_CALLABLE_MEMBER
        Particle(float x,float y, float z, float Ex, float Ey, float Ez, float Z, float amu, float charge)
		{
    
		this->xprevious = x;
		this->yprevious = y;
		this->zprevious = z;
        this->x = x;
        this->y = y;
        this->z = z;
		this->Z = Z;
	    this->charge = charge;
        this->amu = amu;
		this->hitWall = 0.0;
        this->wallIndex = 0;
		this->vx = Ex/std::abs(Ex)*sqrt(2.0*std::abs(Ex)*1.60217662e-19/(amu*1.6737236e-27));
		this->vy = Ey/std::abs(Ey)*sqrt(2.0*std::abs(Ey)*1.60217662e-19/(amu*1.6737236e-27));
		this->vz = Ez/std::abs(Ez)*sqrt(2.0*std::abs(Ez)*1.60217662e-19/(amu*1.6737236e-27));
		    
		if(Ex == 0.0) this->vx = 0.0;
		if(Ey == 0.0) this->vy = 0.0;
		if(Ez == 0.0) this->vz = 0.0;
        };
};
#endif
