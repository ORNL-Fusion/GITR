#ifndef _TESTROUTINE_
#define _TESTROUTINE_

//#ifdef __CUDACC__
//#define CUDA_CALLABLE_MEMBER __host__ __device__
//#define CUDA_CALLABLE_MEMBER_DEVICE __device__
//#else
//#define CUDA_CALLABLE_MEMBER
//#define CUDA_CALLABLE_MEMBER_DEVICE
//#endif

#include <algorithm>
#include "Particle.h"
#include "Boundary.h"
struct test_routine { 
//#ifdef __CUDACC__
//    thrust::device_vector<double>& gridxp;
//#else
   double x;
   double y;
   double z;
   int nx;
   int nz;
   double* gridxp;
   double* gridzp;
   double* datap; 
//#endif

//#ifdef __CUDACC__
//    move_boris(thrust::) : span(_span), boundaryVector(_boundaryVector), nLines(_nLines) {}
//#else
    test_routine(double _x, double _y, double _z,int _nx,int _nz, double* _gridxp,double* _gridzp,
                                double* _datap) : x(_x), y(_y), z(_z),nx(_nx),nz(_nz), gridxp(_gridxp), gridzp(_gridzp),
                                datap(_datap) {}
//#endif    

//CUDA_CALLABLE_MEMBER    
void operator()(double &d) const { 
    std::cout << "testroutine xyz " << x << y << z << std::endl;
    d = interp2dCombined(x,y,z,nx,nz,gridxp,gridzp,datap);
    std::cout << "test routine interp value " << d << std::endl;
}
};

#endif
