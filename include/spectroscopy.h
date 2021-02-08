#ifndef _SPECTROSCOP_
#define _SPECTROSCOP_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <cmath>
#include <vector>
#if USE_CUDA >0
//#if __CUDA_ARCH__ < 600
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
                          //#endif __CUDA_ARCH__ < 600
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
    int nY;
    int nZ;
    float *gridX;
    float *gridY;
    float *gridZ;
    double *bins;
    float dt;

    spec_bin(Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, float *_gridX,float *_gridY,float *_gridZ,
           double * _bins, float _dt) : 
        particlesPointer(_particlesPointer), nBins(_nBins),nX(_nX),nY(_nY), nZ(_nZ), gridX(_gridX),gridY(_gridY),gridZ(_gridZ), bins(_bins),
        dt(_dt) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 
//    int indx_X = 0;
//    int indx_Z = 0;
    float dx = 0.0f;
    float dy = 0.0f;
    float dz = 0.0f;
    float x = particlesPointer->xprevious[indx];
    float y = particlesPointer->yprevious[indx];
    float z = particlesPointer->zprevious[indx];
	    //printf ("z %f \n", z);
#if SPECTROSCOPY > 2
    float dim1 = particlesPointer->xprevious[indx];
#else
  #if USECYLSYMM > 0
    float dim1 = std::sqrt(x*x + y*y);
    #else
    float dim1 = x;
    #endif
#endif

    if ((z > gridZ[0]) && (z < gridZ[nZ-1]))
        {
          if((dim1 > gridX[0]) && (dim1 < gridX[nX-1]))
          {
              dx = gridX[1] - gridX[0];
              dz = gridZ[1] - gridZ[0];
#if SPECTROSCOPY < 3
              int indx_X = std::floor((dim1-gridX[0])/dx);
              int indx_Z = std::floor((z-gridZ[0])/dz);
              int indx_Y = 0;
              int nnYY=1;
	    //printf ("indX %i \n", indx_X);
	    //printf ("indZ %i \n", indx_Z);
#else
              if((y > gridY[0]) && (y < gridY[nY-1]))
              { 
              int indx_X = std::floor((dim1-gridX[0])/dx);
              int indx_Z = std::floor((z-gridZ[0])/dz);
              dy = gridY[1] - gridY[0];
              int indx_Y = std::floor((y-gridY[0])/dy);
              if (indx_Y < 0 || indx_Y >= nY) indx_Y = 0;
              int nnYY = nY;
#endif
              if (indx_X < 0 || indx_X >= nX) indx_X = 0;
              if (indx_Z < 0 || indx_Z >= nZ) indx_Z = 0;
              //std::cout << "gridx0 " << gridX[0] << std::endl;
              //std::cout << "gridz0 " << gridZ[0] << std::endl;
              
              //std::cout << "dx " << dx << std::endl;
              //std::cout << "dz " << dz << std::endl;
              //std::cout << "ind x " << indx_X << "ind z " << indx_Z << std::endl;
              int charge = std::floor(particlesPointer->charge[indx]);
              if(particlesPointer->hitWall[indx]== 0.0)
              {
                  float specWeight = particlesPointer->weight[indx];
		  //printf ("Characters: %f \n", specWeight);
#if USE_CUDA >0
              //for 2d
              /*
              atomicAdd(&bins[nBins*nX*nZ + indx_Z*nX + indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
              if(charge < nBins)
              {
                atomicAdd(&bins[charge*nX*nZ + indx_Z*nX + indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
              }
              */
               //for 3d
	       int index = nBins*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX+ indx_X;
		  //printf ("Index %i \n", index);
              atomicAdd1(&bins[index], specWeight);//0*nX*nZ + indx_Z*nZ + indx_X
              if(charge < nBins)
              {
                atomicAdd1(&bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY + indx_Y*nX+ indx_X], 1.0*specWeight);//0*nX*nZ + indx_Z*nZ + indx_X
              }

#else
              bins[nBins*nX*nnYY*nZ + indx_Z*nX*nnYY  +indx_Y*nX +indx_X] = 
	                          bins[nBins*nX*nnYY*nZ + indx_Z*nX*nnYY+ indx_Y*nX + indx_X] + specWeight;
              if(charge < nBins)
              {
                bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY +indx_Y*nX + indx_X] = 
		                  bins[charge*nX*nnYY*nZ + indx_Z*nX*nnYY+indx_Y*nX + indx_X] + specWeight;
              }
#endif
              }
#if SPECTROSCOPY >2
              }
#endif
              }
        }
    }
};

#endif
