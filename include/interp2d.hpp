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
CUDA_CALLABLE_MEMBER

float interp2d ( float x, float z,int nx, int nz,
    float* gridx,float* gridz,float* data );

CUDA_CALLABLE_MEMBER

float interp2dCombined ( float x, float y, float z,int nx, int nz,
    float* gridx,float* gridz,float* data );
CUDA_CALLABLE_MEMBER

float interp3d ( float x, float y, float z,int nx,int ny, int nz,
    float* gridx,float* gridy, float* gridz,float* data );
CUDA_CALLABLE_MEMBER
void interp3dVector (float* field, float x, float y, float z,int nx,int ny, int nz,
        float* gridx,float* gridy,float* gridz,float* datar, float* dataz, float* datat );
CUDA_CALLABLE_MEMBER
void interp2dVector (float* field, float x, float y, float z,int nx, int nz,
float* gridx,float* gridz,float* datar, float* dataz, float* datat );
CUDA_CALLABLE_MEMBER
void interpFieldAlignedVector (float* field, float x, float y, float z,int nx, int nz,
        float* gridx,float* gridz,float* datar, float* dataz, float* datat,
        int nxB, int nzB, float* gridxB,float* gridzB,float* datarB,float* datazB, float* datatB);
CUDA_CALLABLE_MEMBER
float interp1dUnstructured(float samplePoint,int nx, float max_x, float* data,int &lowInd);
CUDA_CALLABLE_MEMBER
float interp1dUnstructured2(float samplePoint,int nx, float *xdata, float* data);
CUDA_CALLABLE_MEMBER
float interp2dUnstructured(float x,float y,int nx,int ny, float *xgrid,float *ygrid, float* data);
#endif
