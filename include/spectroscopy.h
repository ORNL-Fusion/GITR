#ifndef _SPECTROSCOP_
#define _SPECTROSCOP_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <math.h>
#include <vector>
#if USE_CUDA >0
#if __CUDA_ARCH__ < 600
__device__ double atomicAdd1(double* address, double val)
{
    unsigned long long int* address_as_ull =
                        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
      do {
             assumed = old;
             old = atomicCAS(address_as_ull, assumed,
                            __double_as_longlong(val + 
                                __longlong_as_double(assumed)));
                 // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
                      } while (assumed != old);
                 
                          return __longlong_as_double(old);
                          }
                          #endif
//__device__ double atomicAdd(double* address, double val)
//{
//   unsigned long long int* address_as_ull =
//                                              (unsigned long long int*)address;
//   unsigned long long int old = *address_as_ull, assumed;
//   do {
//         assumed = old;
//         old = atomicCAS(address_as_ull, assumed, 
//               __double_as_longlong(val +                                                 __longlong_as_double(assumed)));
//      } while (assumed != old);
//            return __longlong_as_double(old);
//}
#endif

struct spec_bin { 
    Particles *particlesPointer;
    const int nBins;
    int nX;
    int nZ;
    float *gridX;
    float *gridZ;
    float *bins;
    float dt;

    spec_bin(Particles *_particlesPointer, int _nBins,int _nX, int _nZ, float *_gridX,float *_gridZ,
           float * _bins, float _dt) : 
        particlesPointer(_particlesPointer), nBins(_nBins),nX(_nX), nZ(_nZ), gridX(_gridX),gridZ(_gridZ), bins(_bins),
        dt(_dt) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 
//    int indx_X = 0;
//    int indx_Z = 0;
    float dx = 0.0f;
    float dz = 0.0f;
    float z = particlesPointer->zprevious[indx];
#if USECYLSYMM > 0
    float x = particlesPointer->xprevious[indx];
    float y = particlesPointer->yprevious[indx];
    float dim1 = sqrtf(x*x + y*y);
#else
  float dim1 = particlesPointer->xprevious[indx];
#endif
    if ((z > gridZ[0]) && (z < gridZ[nZ-1]))
        {
          if((dim1 > gridX[0]) && (dim1 < gridX[nX-1]))
          {
              dx = gridX[1] - gridX[0];
              dz = gridZ[1] - gridZ[0];
              int indx_X = (floor((dim1-gridX[0])/dx));
              int indx_Z = floor((z-gridZ[0])/dz + 0.5);
              //std::cout << "gridx0 " << gridX[0] << std::endl;
              //std::cout << "gridz0 " << gridZ[0] << std::endl;
              
              //std::cout << "dx " << dx << std::endl;
              //std::cout << "dz " << dz << std::endl;
              //std::cout << "ind x " << indx_X << "ind z " << indx_Z << std::endl;
              int charge = floor(particlesPointer->charge[indx]);
              if(particlesPointer->hitWall[indx]== 0.0)
              {
#if USE_CUDA >0
              atomicAdd(&bins[nBins*nX*nZ + indx_Z*nX + indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
              if(charge < nBins)
              {
                atomicAdd(&bins[charge*nX*nZ + indx_Z*nX + indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
              }
               //for 3d
              atomicAdd1(&bins[nBins*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX+ indx_X], specWeight);//0*nX*nZ + indx_Z*nZ + indx_X
              if(charge < nBins)
              {
                atomicAdd1(&bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY + indx_Y*nX+ indx_X], 1.0*specWeight);//0*nX*nZ + indx_Z*nZ + indx_X
              }

#else
              bins[nBins*nX*nZ + indx_Z*nX + indx_X] = bins[nBins*nX*nZ + indx_Z*nX + indx_X] + 1.0;
              if(charge < nBins)
              {
                bins[charge*nX*nZ + indx_Z*nX + indx_X] = bins[charge*nX*nZ + indx_Z*nX + indx_X] + 1.0;
              }
#endif
              }
          }
        }
    }
};

#endif
