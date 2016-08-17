#ifndef _TESTROUTINECUDA_
#define _TESTROUTINECUDA_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <algorithm>
#include "Particle.h"
#include "Boundary.h"
struct test_routinecuda { 
   float x;
   float y;
   float z;
   int nx;
   int nz;

   float* gridxp;
   float* gridzp;
   float* datap; 
    
   test_routinecuda(float _x, float _y, float _z,int _nx, int _nz,
           float* _gridxp,float* _gridzp,float* _datap) : 
       x(_x), y(_y), z(_z),nx(_nx), nz(_nz), gridxp(_gridxp),
       gridzp(_gridzp), datap(_datap) {}


CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(float &d) const { 
    //float tmp = gridzp[0];
    //d = tmp;
    d = interp2dcuda(x,y,z,nx,nz,gridxp,gridzp,datap);
}
};

#endif
