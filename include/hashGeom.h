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
     #if USE3DTETGEOM > 0
       float jj = indx/nR;
       int j = floor(jj);
       int i = indx - j*nR;
       int xyzIndx = k*nR*nY + indx;
       float x0 = x[i];
       float y0 = y[j];
       float z0 = z[k];
     #else
       float x0 = x[indx];
       float y0 = 0.0;
       float z0 = z[k];
       int xyzIndx = k*nR + indx;
     #endif
                
     for(int l=0; l<nLines; l++)
     {
        float d1 =((x0 - boundary[l].x1)*(x0 - boundary[l].x1)
                +  (y0 - boundary[l].y1)*(y0 - boundary[l].y1)
                +  (z0 - boundary[l].z1)*(z0 - boundary[l].z1));
        float d2 =((x0 - boundary[l].x2)*(x0 - boundary[l].x2)
                +  (y0 - boundary[l].y2)*(y0 - boundary[l].y2)
                +  (z0 - boundary[l].z2)*(z0 - boundary[l].z2));
          #if USE3DTETGEOM > 0
            float d3 =((x0 - boundary[l].x3)*(x0 - boundary[l].x3)
                    +  (y0 - boundary[l].y3)*(y0 - boundary[l].y3)
                    +  (z0 - boundary[l].z3)*(z0 - boundary[l].z3));
          #endif
          float minOf3 = min(d1,d2);
        //std::cout << "min of two " << minOf3 << std::endl;
          #if USE3DTETGEOM > 0
            minOf3 = min(minOf3,d3);
          #endif
          int minIndClose = n_closeGeomElements;
           for(int m=0; m< n_closeGeomElements; m++)
           {
              if(minDist[xyzIndx*n_closeGeomElements + m] > minOf3)
              {
                  minIndClose = minIndClose-1;
              }
           }

           if(minIndClose < n_closeGeomElements)
           {
        //std::cout << "min INd close " << l << std::endl;
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
           }
      /*     
       if(l == nLines - 1)
       {
           for(int o=0;o<n_closeGeomElements;o++)
              {
                  std::cout << closeGeom[xyzIndx*n_closeGeomElements + o] << " ";
              }
           std::cout << std::endl;
       }
       */
     }
                
   }
};

#endif
