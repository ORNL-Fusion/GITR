#ifndef _HASHGEOM_
#define _HASHGEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include "Boundary.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

#include "interpRateCoeff.hpp"

struct hashGeom {
   int k;
   int nLines; 
   Boundary* boundary;
   float* x;
   float* y;
   float* z;
   int n_closeGeomElements;
   float* minDist;
   int* closeGeom;
   int nR;
   int nY;
   int nZ;


   hashGeom(int _k, int _nLines,
                Boundary* _boundary,
                float* _x,
                float* _y, 
                float* _z, 
                int _n_closeGeomElements, float *_minDist, int *_closeGeom,
                int _nR, int _nY, int _nZ)
               : k(_k), nLines(_nLines),boundary(_boundary), x(_x), y(_y), z(_z), 
               n_closeGeomElements(_n_closeGeomElements), 
               minDist(_minDist), closeGeom(_closeGeom), nR(_nR), nY(_nY), nZ(_nZ) {}
    
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(std::size_t indx) const { 
	              //float kk = indx/(nR*nY);
                  //int k = floor(kk);
                  //float jj = (indx - k*nR*nY)/nR;
                  //int j = floor(jj);
                  //int i = indx - k*nR*nY - j*nR;

                    float jj = indx/nR;
                    int j = floor(jj);
                    int i = indx - j*nR;
                    int xyzIndx = k*nR*nY + indx;
                  float x0 = x[i];
                  float y0 = y[j];
                  float z0 = z[k];
                
                  for(int l=0; l<nLines; l++)
                  {
                      if(boundary[l].Z > 0)
                      {
                      //Get distance
            //          std::cout << "distance calcs " << closeGeomGridr[i] << boundaries[l].x1 << 
              //            closeGeomGridy[j] << boundaries[l].y1 << closeGeomGridz[k] << boundaries[l].z1 << std::endl;
                       float d1 =((x0 - boundary[l].x1)*(x0 - boundary[l].x1)
                               +  (y0 - boundary[l].y1)*(y0 - boundary[l].y1)
                               +  (z0 - boundary[l].z1)*(z0 - boundary[l].z1));
                       float d2 =((x0 - boundary[l].x2)*(x0 - boundary[l].x2)
                               +  (y0 - boundary[l].y2)*(y0 - boundary[l].y2)
                               +  (z0 - boundary[l].z2)*(z0 - boundary[l].z2));
                       float d3 =((x0 - boundary[l].x3)*(x0 - boundary[l].x3)
                               +  (y0 - boundary[l].y3)*(y0 - boundary[l].y3)
                               +  (z0 - boundary[l].z3)*(z0 - boundary[l].z3));
                          //is distance less than min
                //       std::cout << "d123 " << d1 << " " << d2 << " " << d3 << std::endl;
                              float minOf3 = min(d1,d2);
                              minOf3 = min(minOf3,d3);
                          int minIndClose = n_closeGeomElements;
                  //            std::cout << "compare minDist1 to minOf3 " << minOf3 << minIndClose <<std::endl;
                          for(int m=0; m< n_closeGeomElements; m++)
                          {
                    //          std::cout << "minDist[m] " <<minDist1[k*nR_closeGeom*nY_closeGeom*n_closeGeomElements + j*nR_closeGeom*n_closeGeomElements + i*n_closeGeomElements + m] << std::endl;
                             if(minDist[xyzIndx*n_closeGeomElements + m] > minOf3)
                             {
                                 minIndClose = minIndClose-1;
                             }
                          }
                      //    std::cout << "endof m loop " << minIndClose << std::endl;

                          if(minIndClose < n_closeGeomElements)
                          {
                              //%shift numbers down
                              for(int n=n_closeGeomElements-1; n>minIndClose; n--)
                              {
                                   minDist[xyzIndx*n_closeGeomElements + n] = 
                                   minDist[xyzIndx*n_closeGeomElements + n-1];  
                              closeGeom[xyzIndx*n_closeGeomElements+ n] =    
                              closeGeom[xyzIndx*n_closeGeomElements + n-1];
                              }
                              minDist[xyzIndx*n_closeGeomElements + minIndClose] = minOf3;
                              closeGeom[xyzIndx*n_closeGeomElements + minIndClose] = l;
                           //   std::cout << " minof 3 " << minOf3 << std::endl;
                         /*
                              for(int o=0;o<10;o++)
                          {
                              std::cout << "updated min and close " << minDist1[k*nR_closeGeom*nY_closeGeom*n_closeGeomElements + j*nR_closeGeom*n_closeGeomElements + i*n_closeGeomElements + o] <<
                                  " " << closeGeom[k*nR_closeGeom*nY_closeGeom*n_closeGeomElements + j*nR_closeGeom*n_closeGeomElements + i*n_closeGeomElements + o] << std::endl;
                          }
                         */
                          }
                      }
                      }
                
                }
};

#endif
