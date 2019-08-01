#ifndef _BOUNDARY_
#define _BOUNDARY_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <cstdlib>
#include "math.h"
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
template <typename T>
CUDA_CALLABLE_MEMBER
int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

class Boundary 
{
  public:
    int periodic;
    int pointLine;
    int surfaceNumber;
    int surface;
    int inDir;
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float a;
    float b;
    float c;
    float d;
    float plane_norm; //16
    #if USE3DTETGEOM > 0
      float x3;
      float y3;
      float z3;
      float area;
    #else
      float slope_dzdx;
      float intercept_z;
    #endif     
    float Z;
    float amu;
    float potential;
    float ChildLangmuirDist;
    #ifdef __CUDACC__
    //curandState streams[7];
    #else
    //std::mt19937 streams[7];
    #endif
	
    float hitWall;
    float length;
    float distanceToParticle;
    float angle;
    float fd;
    float density;
    float ti;
    float ne;
    float te;
    float debyeLength;
    float larmorRadius;
    float flux;
    float startingParticles;
    float impacts;
    float redeposit;

    CUDA_CALLABLE_MEMBER
    void getSurfaceParallel(float A[],float y,float x)
    {
#if USE3DTETGEOM > 0
        float norm = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        A[1] = (y2-y1)/norm;
#else
        float norm = sqrt((x2-x1)*(x2-x1) + (z2-z1)*(z2-z1));
        A[1] = 0.0;
#endif
        //std::cout << "surf par calc " << x2 << " " << x1 << " " << norm << std::endl;
        A[0] = (x2-x1)/norm;
        A[2] = (z2-z1)/norm;
#if USE3DTETGEOM > 0
#else
#if USECYLSYMM > 0
	    float theta = atan2f(y,x);
	    float B[3] = {0.0f};
            B[0] = cosf(theta)*A[0] - sinf(theta)*A[1];
            B[1] = sinf(theta)*A[0] + cosf(theta)*A[1];
	    A[0] = B[0];
	    A[1] = B[1];
#endif
#endif

    }
    
    CUDA_CALLABLE_MEMBER
    void getSurfaceNormal(float B[],float y,float x)
    {
#if USE3DTETGEOM > 0
        B[0] = -a/plane_norm;
        B[1] = -b/plane_norm;
        B[2] = -c/plane_norm;
#else
        float perpSlope = 0.0;
        if(slope_dzdx == 0.0){perpSlope = 1.0e12;}
        else{
          perpSlope = -sgn(slope_dzdx)/abs(slope_dzdx);}
        float Br = 1.0/sqrt(perpSlope*perpSlope+1.0);
        float Bt = 0.0;
        B[2] = sgn(perpSlope)*sqrt(1-Br*Br);
	#if USECYLSYMM > 0
	    float theta = atan2f(y,x);
            B[0] = cosf(theta)*Br - sinf(theta)*Bt;
            B[1] = sinf(theta)*Br + cosf(theta)*Bt;
            #else
            B[0] = Br;
            B[1] = Bt;
	#endif
        //B[0] = -a/plane_norm;
        //B[1] = -b/plane_norm;
        //B[2] = -c/plane_norm;
        //std::cout << "perp x and z comp " << B[0] << " " << B[2] << std::endl;
#endif
    }
    CUDA_CALLABLE_MEMBER
        void transformToSurface(float C[],float y, float x)
        {
            float X[3] = {0.0f};
            float Y[3] = {0.0f};
            float Z[3] = {0.0f};
            float tmp[3] = {0.0f};
            getSurfaceParallel(X,y,x);
            getSurfaceNormal(Z,y,x);
            Y[0] = Z[1]*X[2] - Z[2]*X[1]; 
            Y[1] = Z[2]*X[0] - Z[0]*X[2]; 
            Y[2] = Z[0]*X[1] - Z[1]*X[0];

            tmp[0] = X[0]*C[0] + Y[0]*C[1] + Z[0]*C[2];
            tmp[1] = X[1]*C[0] + Y[1]*C[1] + Z[1]*C[2];
            tmp[2] = X[2]*C[0] + Y[2]*C[1] + Z[2]*C[2];
            C[0] = tmp[0];
            C[1] = tmp[1];
            C[2] = tmp[2];

        }
//        Boundary(float x1,float y1, float z1, float x2, float y2, float z2,float slope, float intercept, float Z, float amu)
//		{
//    
//		this->x1 = x1;
//		this->y1 = y1;
//		this->z1 = z1;
//        this->x2 = x2;
//        this->y2 = y2;
//        this->z2 = z2;        
//#if USE3DTETGEOM > 0
//#else
//        this->slope_dzdx = slope;
//        this->intercept_z = intercept;
//#endif
//		this->Z = Z;
//		this->amu = amu;
//		this->hitWall = 0.0;
//        array1(amu,0.0);
//        };
};
#endif
