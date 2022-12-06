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

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct hashGeom_sheath {
   //int k;
   int nLines; 
   Boundary* boundary;
   gitr_precision* x;
   gitr_precision* y;
   gitr_precision* z;
   int n_closeGeomElements;
   //gitr_precision* minDist;
   int* closeGeom;
   int nR;
   int nY;
   int nZ;
   int use_3d_geom;


   hashGeom_sheath(int _nLines,
                Boundary* _boundary,
                gitr_precision* _x,
                gitr_precision* _y, 
                gitr_precision* _z, 
                int _n_closeGeomElements, //gitr_precision *_minDist, 
                int *_closeGeom,
                int _nR, int _nY, int _nZ, int use_3d_geom_ )
               : nLines(_nLines),boundary(_boundary), x(_x), y(_y), z(_z), 
               n_closeGeomElements(_n_closeGeomElements), 
               //minDist(_minDist), 
               closeGeom(_closeGeom), nR(_nR), nY(_nY), nZ(_nZ),
               use_3d_geom( use_3d_geom_ ) {}
    
    CUDA_CALLABLE_MEMBER_DEVICE 
    void operator()(std::size_t indx) const { 

       gitr_precision x0;
       gitr_precision y0;
       gitr_precision z0;
       gitr_precision kk;

    if( use_3d_geom > 0 )
    {
       kk = indx/(nR*nY);
       int k = std::floor(kk);
       int jjj = indx - k*nR*nY;
       gitr_precision jj = 1.0*jjj/nR;
       int j = std::floor(jj);
       int i = indx - j*nR - k*(nR*nY);
       x0 = x[i];
       y0 = y[j];
       z0 = z[k];
    }
    else
    {
       kk = indx/(nR);
       int k = std::floor(kk);
       int i = indx - k*(nR);
       x0 = x[i];
       y0 = 0.0;
       z0 = z[k];
    }
       gitr_precision A[3] = {0.0,0.0,0.0};
            gitr_precision B[3] = {0.0,0.0,0.0};
            gitr_precision C[3] = {0.0,0.0,0.0};
            gitr_precision AB[3] = {0.0,0.0,0.0};
            gitr_precision AC[3] = {0.0,0.0,0.0};
            gitr_precision BC[3] = {0.0,0.0,0.0};
            gitr_precision CA[3] = {0.0,0.0,0.0};
            gitr_precision p[3] = {0.0,0.0,0.0};
            gitr_precision Ap[3] = {0.0,0.0,0.0};
            gitr_precision Bp[3] = {0.0,0.0,0.0};
            gitr_precision Cp[3] = {0.0,0.0,0.0};
            gitr_precision normalVector[3] = {0.0,0.0,0.0};
            gitr_precision crossABAp[3] = {0.0,0.0,0.0};
            gitr_precision crossBCBp[3] = {0.0,0.0,0.0};
            gitr_precision crossCACp[3] = {0.0,0.0,0.0};
            gitr_precision signDot0 = 0.0;
            gitr_precision signDot1 = 0.0;
            gitr_precision signDot2 = 0.0;
            gitr_precision totalSigns = 0.0;
#if USE_CUDA
           gitr_precision *minDist  = new gitr_precision[n_closeGeomElements];
           for(int i1=0;i1<n_closeGeomElements;i1++){ minDist[i1] = 1.0e6;}
           //gitr_precision minDist[10] = {1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6,1.0e6};
#else
           sim::Array<gitr_precision> minDist(n_closeGeomElements,1e6);      
#endif
                
    for(int l=0; l<nLines; l++)
    {
        if(boundary[l].Z > 0)
        {
                       gitr_precision a = boundary[l].a;
                       gitr_precision b = boundary[l].b;
                       gitr_precision c = boundary[l].c;
                       gitr_precision d = boundary[l].d;
      gitr_precision perpDist;
    if( use_3d_geom > 0 )
    {
      gitr_precision plane_norm = boundary[l].plane_norm;
      gitr_precision t = -(a*x0 + b*y0 + c*z0 + d)/(a*a + b*b + c*c);
      p[0] = a*t + x0;
      p[1] = b*t + y0;
      p[2] = c*t + z0;
      perpDist = std::sqrt((x0-p[0])*(x0-p[0]) + (y0-p[1])*(y0-p[1]) + (z0-p[2])*(z0-p[2]));
    }
      vectorAssign(boundary[l].x1, boundary[l].y1, 
          boundary[l].z1, A);    
      vectorAssign(boundary[l].x2, boundary[l].y2, 
          boundary[l].z2, B);    
    if( use_3d_geom > 0 )
    {
      vectorAssign(boundary[l].x3, boundary[l].y3, 
          boundary[l].z3, C); 
    }
      vectorSubtract(B,A,AB);
    if( use_3d_geom > 0 )
    {
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
    }
        p[0] = x0;
        p[1] = y0;
        p[2] = z0;
        gitr_precision pA[3] = {0.0};
        gitr_precision cEdge1[3] = {0.0};
        gitr_precision dEdge1[3] = {0.0};
        vectorSubtract(A, p, pA);
        gitr_precision cEdge1mag = vectorDotProduct(pA, AB) / vectorDotProduct(AB, AB);
        gitr_precision distE1 = 1.0e6;
        if (cEdge1mag < 0.0 && cEdge1mag > -1.0) {
          vectorScalarMult(cEdge1mag, AB, cEdge1);
          vectorSubtract(pA, cEdge1, dEdge1);
          distE1 = std::sqrt(vectorDotProduct(dEdge1, dEdge1));
        }
        gitr_precision minEdge;
    if( use_3d_geom > 0 )
    {
        gitr_precision pB[3] = {0.0};
        gitr_precision cEdge2[3] = {0.0};
        gitr_precision dEdge2[3] = {0.0};
        vectorSubtract(B, p, pB);
        gitr_precision cEdge2mag = vectorDotProduct(pB, BC) / vectorDotProduct(BC, BC);
        gitr_precision distE2 = 1.0e6;
        if (cEdge2mag < 0.0 && cEdge2mag > -1.0) {
          vectorScalarMult(cEdge2mag, BC, cEdge2);
          vectorSubtract(pB, cEdge2, dEdge2);
          distE2 = std::sqrt(vectorDotProduct(dEdge2, dEdge2));
        }
        gitr_precision pC[3] = {0.0};
        gitr_precision cEdge3[3] = {0.0};
        gitr_precision dEdge3[3] = {0.0};
        vectorSubtract(C, p, pC);
        gitr_precision cEdge3mag = vectorDotProduct(pC, CA) / vectorDotProduct(CA, CA);
        gitr_precision distE3 = 1.0e6;
        if (cEdge3mag < 0.0 && cEdge3mag > -1.0) {
          vectorScalarMult(cEdge3mag, CA, cEdge3);
          vectorSubtract(pC, cEdge3, dEdge3);
          distE3 = std::sqrt(vectorDotProduct(dEdge3, dEdge3));
        }
        minEdge = std::min(distE1, distE2);
        minEdge = std::min(distE3, minEdge);
    }
    else
    {
          minEdge = distE1;
    }
        gitr_precision d1 =std::sqrt((x0 - boundary[l].x1)*(x0 - boundary[l].x1)
                +  (y0 - boundary[l].y1)*(y0 - boundary[l].y1)
                +  (z0 - boundary[l].z1)*(z0 - boundary[l].z1));
        gitr_precision d2 =std::sqrt((x0 - boundary[l].x2)*(x0 - boundary[l].x2)
                +  (y0 - boundary[l].y2)*(y0 - boundary[l].y2)
                +  (z0 - boundary[l].z2)*(z0 - boundary[l].z2));
            gitr_precision d3;
    if( use_3d_geom > 0 )
    {
            /* Ahoy, Captain! */
            d3 =std::sqrt((x0 - boundary[l].x3)*(x0 - boundary[l].x3)
                    +  (y0 - boundary[l].y3)*(y0 - boundary[l].y3)
                    +  (z0 - boundary[l].z3)*(z0 - boundary[l].z3));
    }
          gitr_precision minOf3 = std::min(d1,d2);
          minOf3 = std::min(minOf3,minEdge);
        //std::cout << "min of two " << minOf3 << std::endl;
    if( use_3d_geom > 0 )
    {
          minOf3 = std::min(minOf3,perpDist);
            minOf3 = std::min(minOf3,d3);
    }
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
