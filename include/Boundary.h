#ifndef _BOUNDARY_
#define _BOUNDARY_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>

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

class Boundary {
	public:
	    float x1;
	    float y1;
	    float z1;
        float x2;
        float y2;
        float z2;
#if USE3DTETGEOM > 0
        float x3;
        float y3;
        float z3;
        float a;
        float b;
        float c;
        float d;
        float plane_norm;
        float area;
#else
        float slope_dzdx;
        float intercept_z;
#endif     
        float Z;
      	float amu;
	#ifdef __CUDACC__
	curandState streams[7];
	#else
        std::mt19937 streams[7];
        #endif
	
	float hitWall;
    float length;
    float distanceToParticle;
    int periodic;
    int pointLine;
    float angle;
    float fd;
    float density;
    float ti;
    float debyeLength;
    float larmorRadius;
    float impacts;
    CUDA_CALLABLE_MEMBER
        Boundary() {
                   };

	CUDA_CALLABLE_MEMBER
        Boundary(float x1,float y1, float z1, float x2, float y2, float z2,float slope, float intercept, float Z, float amu)
		{
    
		this->x1 = x1;
		this->y1 = y1;
		this->z1 = z1;
        this->x2 = x2;
        this->y2 = y2;
        this->z2 = z2;        
#if USE3DTETGEOM > 0
#else
        this->slope_dzdx = slope;
        this->intercept_z = intercept;
#endif
		this->Z = Z;
		this->amu = amu;
		this->hitWall = 0.0;

        };
};
#endif
