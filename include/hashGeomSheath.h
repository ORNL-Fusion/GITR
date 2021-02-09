#ifndef _HASHGEOMSHEATH_
#define _HASHGEOMSHEATH_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include "Boundary.h"
#include "boris.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#include <curand_kernel.h>
#endif

#ifdef __GNUC__ 
#include <random>
#endif

#include "interpRateCoeff.hpp"
#include <cmath>

struct hashGeom_sheath {
   //int k;
   int nLines; 
   Boundary* boundary;
   float* x;
   float* y;
   float* z;
   int n_closeGeomElements;
   //float* minDist;
   int* closeGeom;
   int nR;
   int nY;
   int nZ;


   hashGeom_sheath(int _nLines,
                Boundary* _boundary,
                float* _x,
                float* _y, 
                float* _z, 
                int _n_closeGeomElements, //float *_minDist, 
                int *_closeGeom,
                int _nR, int _nY, int _nZ)
               : nLines(_nLines),boundary(_boundary), x(_x), y(_y), z(_z), 
               n_closeGeomElements(_n_closeGeomElements), 
               //minDist(_minDist), 
               closeGeom(_closeGeom), nR(_nR), nY(_nY), nZ(_nZ) {}
    
    CUDA_CALLABLE_MEMBER_DEVICE 
    void operator()(std::size_t indx) const { 
      #if USE3DTETGEOM > 0
       float kk = indx/(nR*nY);
       int k = std::floor(kk);
       int jjj = indx - k*nR*nY;
       float jj = 1.0*jjj/nR;
       int j = std::floor(jj);
       int i = indx - j*nR - k*(nR*nY);
       int xyzIndx = indx;
       float x0 = x[i];
       float y0 = y[j];
       float z0 = z[k];
      #else
       float kk = indx/(nR);
       int k = std::floor(kk);
       int i = indx - k*(nR);
       float x0 = x[i];
       float y0 = 0.0;
       float z0 = z[k];
       int xyzIndx = indx;
      #endif
       float A[3] = {0.0,0.0,0.0};
            float B[3] = {0.0,0.0,0.0};
            float C[3] = {0.0,0.0,0.0};
            float AB[3] = {0.0,0.0,0.0};
            float AC[3] = {0.0,0.0,0.0};
            float BC[3] = {0.0,0.0,0.0};
            float CA[3] = {0.0,0.0,0.0};
            float p[3] = {0.0,0.0,0.0};
            float Ap[3] = {0.0,0.0,0.0};
            float Bp[3] = {0.0,0.0,0.0};
            float Cp[3] = {0.0,0.0,0.0};
            float normalVector[3] = {0.0,0.0,0.0};
            float crossABAp[3] = {0.0,0.0,0.0};
            float crossBCBp[3] = {0.0,0.0,0.0};
            float crossCACp[3] = {0.0,0.0,0.0};
            float signDot0 = 0.0;
            float signDot1 = 0.0;
            float signDot2 = 0.0;
            float totalSigns = 0.0;
#if USE_CUDA
           float *minDist  = new float[n_closeGeomElements];
           for(int i1=0;i1<n_closeGeomElements;i1++){ minDist[i1] = 1.0e6;}
           //float minDist[10] = {1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6};
#else
           sim::Array<float> minDist(n_closeGeomElements,1e6);      
#endif
                
    for(int l=0; l<nLines; l++)
    {
        if(boundary[l].Z > 0)
        {
                       float a = boundary[l].a;
                       float b = boundary[l].b;
                       float c = boundary[l].c;
                       float d = boundary[l].d;
    #if USE3DTETGEOM > 0
      float plane_norm = boundary[l].plane_norm;
      float t = -(a*x0 + b*y0 + c*z0 + d)/(a*a + b*b + c*c);
      p[0] = a*t + x0;
      p[1] = b*t + y0;
      p[2] = c*t + z0;
      float perpDist = std::sqrt((x0-p[0])*(x0-p[0]) + (y0-p[1])*(y0-p[1]) + (z0-p[2])*(z0-p[2]));
    #endif
      vectorAssign(boundary[l].x1, boundary[l].y1, 
          boundary[l].z1, A);    
      vectorAssign(boundary[l].x2, boundary[l].y2, 
          boundary[l].z2, B);    
    #if USE3DTETGEOM > 0
      vectorAssign(boundary[l].x3, boundary[l].y3, 
          boundary[l].z3, C); 
    #endif
      vectorSubtract(B,A,AB);
    #if USE3DTETGEOM > 0
      vectorSubtract(C,A,AC);
      vectorSubtract(C,B,BC);
      vectorSubtract(A,C,CA);

      vectorSubtract(p,A,Ap);
      vectorSubtract(p,B,Bp);
      vectorSubtract(p,C,Cp);

      vectorCrossProduct(AB,AC,normalVector);
      vectorCrossProduct(AB,Ap,crossABAp);
      vectorCrossProduct(BC,Bp,crossBCBp);
      vectorCrossProduct(CA,Cp,crossCACp);

        signDot0 = std::copysign(1.0,vectorDotProduct(crossABAp, normalVector));
        signDot1 = std::copysign(1.0,vectorDotProduct(crossBCBp, normalVector));
        signDot2 = std::copysign(1.0,vectorDotProduct(crossCACp, normalVector));
        totalSigns = std::abs(signDot0 + signDot1 + signDot2);
        if (totalSigns == 3.0) {
        } else
          perpDist = 1.0e6;
#endif
        p[0] = x0;
        p[1] = y0;
        p[2] = z0;
        float pA[3] = {0.0};
        float cEdge1[3] = {0.0};
        float dEdge1[3] = {0.0};
        vectorSubtract(A, p, pA);
        float cEdge1mag = vectorDotProduct(pA, AB) / vectorDotProduct(AB, AB);
        float distE1 = 1.0e6;
        if (cEdge1mag < 0.0 && cEdge1mag > -1.0) {
          vectorScalarMult(cEdge1mag, AB, cEdge1);
          vectorSubtract(pA, cEdge1, dEdge1);
          distE1 = std::sqrt(vectorDotProduct(dEdge1, dEdge1));
        }
#if USE3DTETGEOM > 0
        float pB[3] = {0.0};
        float cEdge2[3] = {0.0};
        float dEdge2[3] = {0.0};
        vectorSubtract(B, p, pB);
        float cEdge2mag = vectorDotProduct(pB, BC) / vectorDotProduct(BC, BC);
        float distE2 = 1.0e6;
        if (cEdge2mag < 0.0 && cEdge2mag > -1.0) {
          vectorScalarMult(cEdge2mag, BC, cEdge2);
          vectorSubtract(pB, cEdge2, dEdge2);
          distE2 = std::sqrt(vectorDotProduct(dEdge2, dEdge2));
        }
        float pC[3] = {0.0};
        float cEdge3[3] = {0.0};
        float dEdge3[3] = {0.0};
        vectorSubtract(C, p, pC);
        float cEdge3mag = vectorDotProduct(pC, CA) / vectorDotProduct(CA, CA);
        float distE3 = 1.0e6;
        if (cEdge3mag < 0.0 && cEdge3mag > -1.0) {
          vectorScalarMult(cEdge3mag, CA, cEdge3);
          vectorSubtract(pC, cEdge3, dEdge3);
          distE3 = std::sqrt(vectorDotProduct(dEdge3, dEdge3));
        }
        float minEdge = std::min(distE1, distE2);
        minEdge = std::min(distE3, minEdge);
#else
          //
          float minEdge = distE1;
#endif
        float d1 =std::sqrt((x0 - boundary[l].x1)*(x0 - boundary[l].x1)
                +  (y0 - boundary[l].y1)*(y0 - boundary[l].y1)
                +  (z0 - boundary[l].z1)*(z0 - boundary[l].z1));
        float d2 =std::sqrt((x0 - boundary[l].x2)*(x0 - boundary[l].x2)
                +  (y0 - boundary[l].y2)*(y0 - boundary[l].y2)
                +  (z0 - boundary[l].z2)*(z0 - boundary[l].z2));
          #if USE3DTETGEOM > 0
            float d3 =std::sqrt((x0 - boundary[l].x3)*(x0 - boundary[l].x3)
                    +  (y0 - boundary[l].y3)*(y0 - boundary[l].y3)
                    +  (z0 - boundary[l].z3)*(z0 - boundary[l].z3));
          #endif
          float minOf3 = std::min(d1,d2);
          minOf3 = std::min(minOf3,minEdge);
        //std::cout << "min of two " << minOf3 << std::endl;
          #if USE3DTETGEOM > 0
          minOf3 = std::min(minOf3,perpDist);
            minOf3 = std::min(minOf3,d3);
          #endif
          int minIndClose = n_closeGeomElements;
           for(int m=0; m< n_closeGeomElements; m++)
           {
              if(minDist[m] > minOf3)
              {
                  minIndClose = minIndClose-1;
              }
           }

           if((minIndClose < n_closeGeomElements) && (minIndClose > -1))
           {
               //%shift numbers down
               for(int n=n_closeGeomElements-1; n>minIndClose; n--)
               {
                    minDist[n] = 
                    minDist[n-1];  
               closeGeom[indx*n_closeGeomElements+ n] =    
               closeGeom[indx*n_closeGeomElements + n-1];
               }
               minDist[minIndClose] = minOf3;
              closeGeom[indx*n_closeGeomElements + minIndClose] = l;
     }
        }
    }
#if USE_CUDA
     delete[] minDist;
#endif
                
                }
};

#endif
