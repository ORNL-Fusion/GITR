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
#include "array.h"
#include "managed_allocation.h"

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
        float potential;
        float ChildLangmuirDist;
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
    float startingParticles;
    float impacts;
    float redeposit;
    int nE = 1000;
    int nA = 90;
    float E0 = 0.0;
    float E1 = 1000.0;
    float A0 = 0.0;
    float A1 = 90.0;
    //sim::Array<float> array2;
    std::vector<float> array1;// = std::vector<float>(90000);
    float array3[90000] = {0.0};
    CUDA_CALLABLE_MEMBER
        Boundary() : array1(1000,0.0) {
                   };
//  CUDA_CALLABLE_MEMBER
//    Boundary(std::size_t nP,int nE) : array1(nE,0.0,1) {};

  //void initArray(int nE)
  //{
  //    this->array1(nE,0.0);
  //}
  Boundary & operator=(const Boundary &rhs) 
  {   //array1(6);
      //array2 = rhs.array2;
//      //x1 = rhs.x1;
        array1 = rhs.array1;
        //array2 = rhs.array2;
//          std::cout << "here 1" << std::endl;
       //array2.resize(rhs.array2.size());
//     // for(int i=0;i<rhs.array1.size();i++)
//     // {
//     //     std::cout << "here " << std::endl;
//     //    array1[i] = rhs.array1[i];
//     // }
//      //array1(10);
      //std::swap(array2,rhs.array2);
      return *this;
  };
  //CUDA_CALLABLE_MEMBER
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
