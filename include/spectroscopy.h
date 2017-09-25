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

struct spec_bin { 
    Particles *particlesPointer;
    const int nBins;
    int nX;
    int nY;
    int nZ;
    float *gridX;
    float *gridY;
    float *gridZ;
    float *bins;
    float dt;

    spec_bin(Particles *_particlesPointer, int _nBins,int _nX,int _nY, int _nZ, float *_gridX,float *_gridY,float *_gridZ,
           float * _bins, float _dt) : 
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
#if USECYLSYMM > 0
    //float dim1 = sqrtf(x*x + y*y);
    float dim1 = particlesPointer->xprevious[indx];
#else
  float dim1 = particlesPointer->xprevious[indx];
#endif
    if ((z > gridZ[0]) && (z < gridZ[nZ-1]))
        {
          if((dim1 > gridX[0]) && (dim1 < gridX[nX-1]))
          {
              dx = gridX[1] - gridX[0];
              dy = gridY[1] - gridY[0];
              dz = gridZ[1] - gridZ[0];
#if USECYLSYMM > 0
              int indx_X = floor((dim1-gridX[0])/dx);
              int indx_Z = floor((z-gridZ[0])/dz);
#else
              int indx_X = floor((dim1-gridX[0])/dx+0.5);
              int indx_Z = floor((z-gridZ[0])/dz + 0.5);
#endif
              int indx_Y = floor((y-gridY[0])/dy);
              if (indx_X < 0 || indx_X >= nX) indx_X = 0;
              if (indx_Y < 0 || indx_Y >= nY) indx_Y = 0;
              if (indx_Z < 0 || indx_Z >= nZ) indx_Z = 0;
              //std::cout << "gridx0 " << gridX[0] << std::endl;
              //std::cout << "gridz0 " << gridZ[0] << std::endl;
              
              //std::cout << "dx " << dx << std::endl;
              //std::cout << "dz " << dz << std::endl;
              //std::cout << "ind x " << indx_X << "ind z " << indx_Z << std::endl;
              int charge = floor(particlesPointer->charge[indx]);
              if(particlesPointer->hitWall[indx]== 0.0)
              {
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
              atomicAdd(&bins[nBins*nX*nY*nZ + indx_Z*nX*nY +indx_Y*nX+ indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
              if(charge < nBins)
              {
                atomicAdd(&bins[charge*nX*nY*nZ + indx_Z*nX*nY + indx_Y*nX+ indx_X], 1.0);//0*nX*nZ + indx_Z*nZ + indx_X
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
