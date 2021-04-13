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
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
struct test_routinecuda { 
   gitr_precision x;
   gitr_precision y;
   gitr_precision z;
   int nx;
   int nz;

   gitr_precision* gridxp;
   gitr_precision* gridzp;
   gitr_precision* datap; 
    
   test_routinecuda(gitr_precision _x, gitr_precision _y, gitr_precision _z,int _nx, int _nz,
           gitr_precision* _gridxp,gitr_precision* _gridzp,gitr_precision* _datap) : 
       x(_x), y(_y), z(_z),nx(_nx), nz(_nz), gridxp(_gridxp),
       gridzp(_gridzp), datap(_datap) {}


CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(gitr_precision &d) const { 
    //gitr_precision tmp = gridzp[0];
    //d = tmp;
    d = interp2dCombined(x,y,z,nx,nz,gridxp,gridzp,datap);
}
};

#endif
