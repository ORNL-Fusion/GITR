#ifndef _INTERP2D_
#define _INTERP2D_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <thrust/device_vector.h>
#include <vector>
#include "math.h"
//using namespace std;

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
CUDA_CALLABLE_MEMBER

gitr_precision interp2d ( gitr_precision x, gitr_precision z,int nx, int nz,
    gitr_precision* gridx,gitr_precision* gridz,gitr_precision* data );

CUDA_CALLABLE_MEMBER

gitr_precision interp2dCombined ( gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz,
    gitr_precision* gridx,gitr_precision* gridz,gitr_precision* data);
CUDA_CALLABLE_MEMBER

gitr_precision interp3d ( gitr_precision x, gitr_precision y, gitr_precision z,int nx,int ny, int nz,
    gitr_precision* gridx,gitr_precision* gridy, gitr_precision* gridz,gitr_precision* data );
CUDA_CALLABLE_MEMBER

gitr_precision interp3d_nearest ( gitr_precision x, gitr_precision y, gitr_precision z,int nx,int ny, int nz,
    gitr_precision* gridx,gitr_precision* gridy, gitr_precision* gridz,gitr_precision* data );
CUDA_CALLABLE_MEMBER
void interp3dVector (gitr_precision* field, gitr_precision x, gitr_precision y, gitr_precision z,int nx,int ny, int nz,
        gitr_precision* gridx,gitr_precision* gridy,gitr_precision* gridz,gitr_precision* datar, gitr_precision* dataz, gitr_precision* datat );
CUDA_CALLABLE_MEMBER
void interp2dVector (gitr_precision* field, gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz,
gitr_precision* gridx,gitr_precision* gridz,gitr_precision* datar, gitr_precision* dataz, gitr_precision* datat);
CUDA_CALLABLE_MEMBER
void interpFieldAlignedVector ( gitr_precision* field, gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz,
        gitr_precision* gridx,gitr_precision* gridz,gitr_precision* datar, gitr_precision* dataz, gitr_precision* datat,
        int nxB, int nzB, gitr_precision* gridxB,gitr_precision* gridzB,gitr_precision* datarB,gitr_precision* datazB, gitr_precision* datatB);
CUDA_CALLABLE_MEMBER
gitr_precision interp1dUnstructured(gitr_precision samplePoint,int nx, gitr_precision max_x, gitr_precision* data,int &lowInd);
CUDA_CALLABLE_MEMBER
gitr_precision interp1dUnstructured2(gitr_precision samplePoint,int nx, gitr_precision *xdata, gitr_precision* data);
CUDA_CALLABLE_MEMBER
gitr_precision interp2dUnstructured(gitr_precision x,gitr_precision y,int nx,int ny, gitr_precision *xgrid,gitr_precision *ygrid, gitr_precision* data);
#endif
