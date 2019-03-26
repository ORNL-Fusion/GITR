#ifndef _BORIS_
#define _BORIS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#include "thrust/extrema.h"
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <algorithm>
#include "Particles.h"
#include "Boundary.h"
#include "interp2d.hpp"
#if USE_BOOST
#include <boost/timer/timer.hpp>
using namespace boost::timer;
#endif
//template <typename T>
//CUDA_CALLABLE_MEMBER
//int sgn(T val) {
//            return (T(0) < val) - (val < T(0));
//}

CUDA_CALLABLE_MEMBER
void vectorAdd(float A[], float B[],float C[])
{
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
}

CUDA_CALLABLE_MEMBER
void vectorSubtract(float A[], float B[],float C[])
{
    C[0] = A[0] - B[0];
    C[1] = A[1] - B[1];
    C[2] = A[2] - B[2];
}

CUDA_CALLABLE_MEMBER
void vectorScalarMult(float a, float B[],float C[])
{
    C[0] = a*B[0];
    C[1] = a*B[1];
    C[2] = a*B[2];
}

CUDA_CALLABLE_MEMBER
void vectorAssign(float a, float b,float c, float D[])
{
    D[0] = a;
    D[1] = b;
    D[2] = c;
}

CUDA_CALLABLE_MEMBER
float vectorNorm(float A[])
{
    float norm = 0.0f;
    norm = sqrtf(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);

        return norm;
}
CUDA_CALLABLE_MEMBER
void vectorNormalize(float A[],float B[])
{
    float norm = 0.0f;
    norm = sqrtf(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    B[0] = A[0]/norm;
    B[1] = A[1]/norm;
    B[2] = A[2]/norm;

}

CUDA_CALLABLE_MEMBER
float vectorDotProduct(float A[], float B[])
{
    float c = A[0]*B[0] +  A[1]*B[1] + A[2]*B[2];
    return c;
}

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(float A[], float B[], float C[])
{
    float tmp[3] = {0.0f,0.0f,0.0f};
    tmp[0] = A[1]*B[2] - A[2]*B[1];
    tmp[1] = A[2]*B[0] - A[0]*B[2];
    tmp[2] = A[0]*B[1] - A[1]*B[0];

    C[0] = tmp[0];
    C[1] = tmp[1];
    C[2] = tmp[2];
}
CUDA_CALLABLE_MEMBER

float getE ( float x0, float y, float z, float E[], Boundary *boundaryVector, int nLines,
       int nR_closeGeom, int nY_closeGeom,int nZ_closeGeom, int n_closeGeomElements, 
       float *closeGeomGridr,float *closeGeomGridy, float *closeGeomGridz, int *closeGeom, 
         int&  closestBoundaryIndex) {
#if USE3DTETGEOM > 0
    float Emag = 0.0f;
    float Er = 0.0f;
    float Et = 0.0f;
      float p0[3] = {x0,y,z};
    float angle = 0.0f;
	float fd = 0.0f;
	float pot = 0.0f;
      float a = 0.0;
      float b = 0.0;
      float c = 0.0;
      float d = 0.0;
      float plane_norm = 0.0;
      float pointToPlaneDistance0 = 0.0;
      float pointToPlaneDistance1 = 0.0;
      float signPoint0 = 0.0;
      float signPoint1 = 0.0;
      float t = 0.0;
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
      float p0A[3] = {0.0,0.0,0.0};
      float p0B[3] = {0.0,0.0,0.0};
      float p0C[3] = {0.0,0.0,0.0};
      float p0AB[3] = {0.0,0.0,0.0};
      float p0BC[3] = {0.0,0.0,0.0};
      float p0CA[3] = {0.0,0.0,0.0};
      float p0Anorm = 0.0f;
      float p0Bnorm = 0.0f;
      float p0Cnorm = 0.0f;
      float normalVector[3] = {0.0,0.0,0.0};
      float crossABAp[3] = {0.0,0.0,0.0};
      float crossBCBp[3] = {0.0,0.0,0.0};
      float crossCACp[3] = {0.0,0.0,0.0};
      float directionUnitVector[3] = {0.0f,0.0f,0.0f};
      float dot0 = 0.0f;
      float dot1 = 0.0f;
      float dot2 = 0.0f;

      float normAB = 0.0f;
      float normBC = 0.0f;
      float normCA = 0.0f;
      float ABhat[3] = {0.0f,0.0f,0.0f};
      float BChat[3] = {0.0f,0.0f,0.0f};
      float CAhat[3] = {0.0f,0.0f,0.0f};
      float tAB = 0.0f;
      float tBC = 0.0f;
      float tCA = 0.0f;
      float projP0AB[3] = {0.0f,0.0f,0.0f};
      float projP0BC[3] = {0.0f,0.0f,0.0f};
      float projP0CA[3] = {0.0f,0.0f,0.0f};
      float p0ABdist = 0.0f;
      float p0BCdist = 0.0f;
      float p0CAdist = 0.0f;
      float perpDist = 0.0f;
      float signDot0 = 0.0;
      float signDot1 = 0.0;
      float signDot2 = 0.0;
      float totalSigns = 0.0;
      float minDistance = 1e12f;
      int nBoundariesCrossed = 0;
      int boundariesCrossed[6] = {0,0,0,0,0,0};
        int minIndex=0;
      float distances[7] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
      float normals[21] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                           0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                           0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
#if GEOM_HASH_SHEATH > 0
  float dr = closeGeomGridr[1] - closeGeomGridr[0];
  float dy = closeGeomGridy[1] - closeGeomGridy[0];
  float dz = closeGeomGridz[1] - closeGeomGridz[0];
  int rInd = floor((x0 - closeGeomGridr[0])/dr + 0.5f);
  int yInd = floor((y - closeGeomGridy[0])/dy + 0.5f);
  int zInd = floor((z - closeGeomGridz[0])/dz + 0.5f);
  int i;
  if(rInd < 0 || rInd >= nR_closeGeom)
    rInd =0;
  if(yInd < 0 || yInd >= nY_closeGeom)
    yInd =0;
  if(zInd < 0 || zInd >= nZ_closeGeom)
    zInd =0;

  for (int k=0; k< n_closeGeomElements; k++) //n_closeGeomElements
    {
       i = closeGeom[zInd*nY_closeGeom*nR_closeGeom*n_closeGeomElements 
                   + yInd*nR_closeGeom*n_closeGeomElements
                   + rInd*n_closeGeomElements + k];
       //closestBoundaryIndex = i;
       //std::cout << "closest boundaries to check " << i << std::endl;
#else
      for (int i=0; i<nLines; i++)
      {
#endif
            //std::cout << "Z and index " << boundaryVector[i].Z << " " << i << std::endl;
        //if (boundaryVector[i].Z != 0.0)
        //{
            //std::cout << "Z and index " << boundaryVector[i].Z << " " << i << std::endl;
        a = boundaryVector[i].a;
        b = boundaryVector[i].b;
        c = boundaryVector[i].c;
        d = boundaryVector[i].d;
        plane_norm = boundaryVector[i].plane_norm;
        pointToPlaneDistance0 = (a*p0[0] + b*p0[1] + c*p0[2] + d)/plane_norm;
        //std::cout << "abcd plane_norm "<< a  << " " << b << " " << c << " " << d << " " << plane_norm << std::endl;
        //std::cout << i << std::endl;// " point to plane dist "  << pointToPlaneDistance0 << std::endl;
        //pointToPlaneDistance1 = (a*p1[0] + b*p1[1] + c*p1[2] + d)/plane_norm;    
        //signPoint0 = sgn(pointToPlaneDistance0);
        //signPoint1 = sgn(pointToPlaneDistance1);
        vectorAssign(a/plane_norm,b/plane_norm,c/plane_norm,normalVector);
                //vectorNormalize(normalVector,normalVector);
                //std::cout << "normal " << normalVector[0] << " " << normalVector[1] << " " << normalVector[2] << std::endl;
          vectorAssign(p0[0] - pointToPlaneDistance0*normalVector[0],
                       p0[1] - pointToPlaneDistance0*normalVector[1],
                       p0[2] - pointToPlaneDistance0*normalVector[2], p);
        
         vectorAssign(boundaryVector[i].x1, boundaryVector[i].y1,
                     boundaryVector[i].z1, A);
         vectorAssign(boundaryVector[i].x2, boundaryVector[i].y2,
                      boundaryVector[i].z2, B);
         vectorAssign(boundaryVector[i].x3, boundaryVector[i].y3,
                     boundaryVector[i].z3, C);

         vectorSubtract(B,A,AB);
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
       /*  
         dot0 = vectorDotProduct(crossABAp,normalVector);
            dot1 = vectorDotProduct(crossBCBp,normalVector);
            dot2 = vectorDotProduct(crossCACp,normalVector);
         */
            signDot0 = sgn(vectorDotProduct(crossABAp,normalVector));
         signDot1 = sgn(vectorDotProduct(crossBCBp,normalVector));
         signDot2 = sgn(vectorDotProduct(crossCACp,normalVector));
        /*  
         if(dot0 == 0.0) signDot0 = 1;
         if(dot1 == 0.0) signDot1 = 1;
         if(dot2 == 0.0) signDot2 = 1;
         
         if(vectorNorm(crossABAp) == 0.0) signDot0 = 1;
         if(vectorNorm(crossBCBp) == 0.0) signDot1 = 1;
         if(vectorNorm(crossCACp) == 0.0) signDot2 = 1;
         */
         totalSigns = abs(signDot0 + signDot1 + signDot2);
         /*
         if(vectorNorm(crossABAp) == 0.0) totalSigns = 3;
         if(vectorNorm(crossBCBp) == 0.0) totalSigns = 3;
         if(vectorNorm(crossCACp) == 0.0) totalSigns = 3;
         */
         //std::cout << dot0 << " signDot0 " << signDot0 << std::endl;
         //std::cout << dot1 << " signDot1 " << signDot1 << std::endl;
         //std::cout << dot2 << " signDot2 " << signDot2 << std::endl;
         vectorSubtract(A,p0,p0A);
         vectorSubtract(B,p0,p0B);
         vectorSubtract(C,p0,p0C);
         
         p0Anorm = vectorNorm(p0A);   
         p0Bnorm = vectorNorm(p0B);   
         p0Cnorm = vectorNorm(p0C);
         distances[1] = p0Anorm;   
         distances[2] = p0Bnorm;   
         distances[3] = p0Cnorm;   
             normals[3] = p0A[0]/p0Anorm;
             normals[4] = p0A[1]/p0Anorm;
             normals[5] = p0A[2]/p0Anorm;
             normals[6] = p0B[0]/p0Bnorm;
             normals[7] = p0B[1]/p0Bnorm;
             normals[8] = p0B[2]/p0Bnorm;
             normals[9] = p0C[0]/p0Cnorm;
             normals[10] = p0C[1]/p0Cnorm;
             normals[11] = p0C[2]/p0Cnorm;
         //std::cout << "point to plane " << pointToPlaneDistance0 << std::endl;
         //std::cout << "point to ABC " << p0Anorm << " " << p0Bnorm << " " << p0Cnorm << std::endl;
         //std::cout << "total Signs " << totalSigns << std::endl;
         normAB = vectorNorm(AB);
         normBC = vectorNorm(BC);
         normCA = vectorNorm(CA);
         vectorAssign(AB[0]/normAB,AB[1]/normAB,AB[2]/normAB,ABhat);
         vectorAssign(BC[0]/normBC,BC[1]/normBC,BC[2]/normBC,BChat);
         vectorAssign(CA[0]/normCA,CA[1]/normCA,CA[2]/normCA,CAhat);
         
         tAB = vectorDotProduct(p0A,ABhat);
         tBC = vectorDotProduct(p0B,BChat);
         tCA = vectorDotProduct(p0C,CAhat);
         tAB = -1*tAB;
         tBC = -1*tBC;
         tCA = -1*tCA;
         /*
         std::cout << "A " << A[0] << " " << A[1] << " " << A[2] << std::endl;   
         std::cout << "B " << B[0] << " " << B[1] << " " << B[2] << std::endl;   
         std::cout << "C " << C[0] << " " << C[1] << " " << C[2] << std::endl;   
         std::cout << "ABhat " << ABhat[0] << " " << ABhat[1] << " " << ABhat[2] << std::endl; 
         std::cout << "CAhat " << CAhat[0] << " " << CAhat[1] << " " << CAhat[2] << std::endl; 
         std::cout << "p0C " << p0C[0] << " " << p0C[1] << " " << p0C[2] << std::endl; 
         std::cout << "tAB and norm AB " << tAB << " "<<  normAB << std::endl;  
         std::cout << "tBC and norm BC " << tBC << " "<<  normBC << std::endl;  
         std::cout << "tCA and norm CA " << tCA << " "<<  normCA << std::endl;  
         */
         if((tAB > 0.0) && (tAB < normAB))
         {
             vectorScalarMult(tAB,ABhat,projP0AB);
             vectorAdd(A,projP0AB,projP0AB);
             vectorSubtract(projP0AB,p0,p0AB);
             p0ABdist = vectorNorm(p0AB);
             distances[4] = p0ABdist;   
             normals[12] = p0AB[0]/p0ABdist;
             normals[13] = p0AB[1]/p0ABdist;
             normals[14] = p0AB[2]/p0ABdist;

         }
         else
         {
             p0ABdist = 1e12f;
             distances[4] = p0ABdist;   
         } 
         
         
         if((tBC > 0.0) && (tBC < normBC))
         {
             vectorScalarMult(tBC,ABhat,projP0BC);
             vectorAdd(B,projP0BC,projP0BC);
             vectorSubtract(projP0BC,p0,p0BC);
             p0BCdist = vectorNorm(p0BC);
             distances[5] = p0BCdist;   
             normals[15] = p0BC[0]/p0BCdist;
             normals[16] = p0BC[1]/p0BCdist;
             normals[17] = p0BC[2]/p0BCdist;

         }
         else
         {
             p0BCdist = 1e12f;
             distances[5] = p0BCdist;   

         } 
         
         if((tCA > 0.0) && (tCA < normCA))
         {
             vectorScalarMult(tCA,CAhat,projP0CA);
             vectorAdd(C,projP0CA,projP0CA);
             //std::cout << "projP0CA " << projP0CA[0] << " " << projP0CA[1] << " " << projP0CA[2] << std::endl; 
             vectorSubtract(projP0CA,p0,p0CA);
             p0CAdist = vectorNorm(p0CA);
             distances[6] = p0CAdist;   
             normals[18] = p0CA[0]/p0CAdist;
             normals[19] = p0CA[1]/p0CAdist;
             normals[20] = p0CA[2]/p0CAdist;
             //std::cout << "p0CA " << p0CA[0] << " " << p0CA[1] << " " << p0CA[2] << std::endl; 
         }
         else
         {
             p0CAdist = 1e12f;
             distances[6] = p0CAdist;   
         } 

         if (totalSigns == 3.0)
         {
             //if (fabs(pointToPlaneDistance0) < minDistance)
             //{
                perpDist = fabs(pointToPlaneDistance0); 
                //minDistance = fabs(pointToPlaneDistance0);
                //std::cout << "p " << p[0] << " " << p[1] << " " << p[2] << std::endl;
                //std::cout << "p0 " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
                vectorSubtract(p,p0 ,normalVector);
                //std::cout << "unit vec " << directionUnitVector[0] << " " << directionUnitVector[1] << 
                //    " " << directionUnitVector[2] << std::endl;
                vectorNormalize(normalVector,normalVector);
                //std::cout << "unit vec " << directionUnitVector[0] << " " << directionUnitVector[1] << 
                //    " " << directionUnitVector[2] << std::endl;
                //std::cout << "perp distance " << std::endl;
             distances[0] = perpDist;   
             normals[0] = normalVector[0];
             normals[1] = normalVector[1];
             normals[2] = normalVector[2];
             //}
         }
         else
         {
             perpDist = 1e12f;
             distances[0] = perpDist;   
             /*
             if (p0Anorm < p0Bnorm)
             {
                     if (p0Anorm < p0Cnorm)
                     {
                        if (p0Anorm < minDistance)
                        {
                         minDistance = p0Anorm;
                         vectorAssign(p0A[0]/p0Anorm,p0A[1]/p0Anorm,p0A[2]/p0Anorm,
                                      directionUnitVector);
                         minIndex = i;
                        std::cout << "p0A " << std::endl;
                        }
                     }
                     else
                     {
                        if (p0Cnorm < minDistance)
                        {
                         minDistance = p0Cnorm;
                         vectorAssign(p0C[0]/p0Cnorm,p0C[1]/p0Cnorm,p0C[2]/p0Cnorm,
                                      directionUnitVector);
                         minIndex = i;
                        std::cout << "p0C " << p0Cnorm << std::endl;
                        }
                     }
             }
             else
             {
                     if (p0Bnorm < p0Cnorm)
                     {
                        if (p0Bnorm < minDistance)
                        {
                         minDistance = p0Bnorm;
                         vectorAssign(p0B[0]/p0Bnorm,p0B[1]/p0Bnorm,p0B[2]/p0Bnorm,
                                      directionUnitVector);
                         minIndex = i;
                        std::cout << "p0B " << std::endl;
                        }
                     }
                     else
                     {
                        if (p0Cnorm < minDistance)
                        {
                         minDistance = p0Cnorm;
                         vectorAssign(p0C[0]/p0Cnorm,p0C[1]/p0Cnorm,p0C[2]/p0Cnorm,
                                      directionUnitVector);
                         minIndex = i;
                        std::cout << "p0C two " << std::endl;
                        }
                     }

             }      
             */
         }
         int index = 0;
         for(int j = 0; j < 7; j++)
         {
            if(distances[j] < distances[index])
            index = j;              
         }

         if (distances[index] < minDistance)
         {
                 minDistance = distances[index];
                 vectorAssign(normals[index*3], normals[index*3+1],normals[index*3+2], directionUnitVector);
                 //std::cout << "min dist " << minDistance << std::endl;
                 //std::cout << "min normal " << normals[index*3] << " " 
                 //   <<normals[index*3+1] << " " << normals[index*3+2] << std::endl;
               //closestBoundaryIndex = i;
          closestBoundaryIndex = i;
          minIndex = i;
         }
         //std::cout << "perp dist " << perpDist << std::endl;
         //std::cout << "point to AB BC CA " << p0ABdist << " " << p0BCdist << " " << p0CAdist << std::endl;
        //}
       }
      //vectorScalarMult(-1.0,directionUnitVector,directionUnitVector);
      //std::cout << "min dist " << minDistance << std::endl;
#else      
                
    float Emag = 0.0f;
	float fd = 0.0f;
	float pot = 0.0f;
    int minIndex = 0;
    float minDistance = 1e12f;
    int direction_type;
    float tol = 1e12f;
    float point1_dist;
    float point2_dist;
    float perp_dist;
    float directionUnitVector[3] = {0.0f,0.0f,0.0f};
    float vectorMagnitude;
    float max = 0.0f;
    float min = 0.0f;
    float angle = 0.0f;
    float Er = 0.0f;
    float Et = 0.0f;
    float Bfabsfperp = 0.0f;
    float distanceToParticle = 0.0f;
    int pointLine=0;
//#if EFIELD_INTERP ==1
#if USECYLSYMM > 0
    float x = sqrtf(x0*x0 + y*y);
#else
    float x = x0;
#endif 

#if GEOM_HASH_SHEATH > 0
  float dr = closeGeomGridr[1] - closeGeomGridr[0];
  float dz = closeGeomGridz[1] - closeGeomGridz[0];
  int rInd = floor((x - closeGeomGridr[0])/dr + 0.5f);
  int zInd = floor((z - closeGeomGridz[0])/dz + 0.5f);
  int j;
  for (int k=0; k< n_closeGeomElements; k++) //n_closeGeomElements
    {
       j = closeGeom[zInd*nR_closeGeom*n_closeGeomElements + rInd*n_closeGeomElements + k];

#else
    for (int j=0; j< nLines; j++)
    {
#endif
        if (boundaryVector[j].Z != 0.0)
        {
            point1_dist = sqrtf((x - boundaryVector[j].x1)*(x - boundaryVector[j].x1) + 
                    (z - boundaryVector[j].z1)*(z - boundaryVector[j].z1));
            point2_dist = sqrtf((x - boundaryVector[j].x2)*(x - boundaryVector[j].x2) + 
                                        (z - boundaryVector[j].z2)*(z - boundaryVector[j].z2));
            perp_dist = (boundaryVector[j].slope_dzdx*x - z + boundaryVector[j].intercept_z)/
                sqrtf(boundaryVector[j].slope_dzdx*boundaryVector[j].slope_dzdx + 1.0f);   

            if (point1_dist > point2_dist)
            {
                max = point1_dist;
                min = point2_dist;
            }
            else
            {
                max = point2_dist;
                min = point1_dist;
            }
    //        std::cout << "p1dist p2dist perpDist " << point1_dist << " " << point2_dist << " " << perp_dist << std::endl;
            if (boundaryVector[j].length*boundaryVector[j].length + perp_dist*perp_dist >=
                    max*max)
            {
                //boundaryVector[j].distanceToParticle =fabsf( perp_dist);
                distanceToParticle = fabsf(perp_dist);
                //boundaryVector[j].pointLine = 1;
                pointLine = 1;
            }
            else
            {
                //boundaryVector[j].distanceToParticle = min;
                distanceToParticle = min;
                if (boundaryVector[j].distanceToParticle == point1_dist)
                {
                    pointLine = 2;
                }
                else
                {
                    pointLine = 3;
                }
            }

            if (distanceToParticle < minDistance)
            {
                minDistance = distanceToParticle;
                minIndex = j;
                closestBoundaryIndex = j;
                direction_type = pointLine;
            }
        }
        else
        {
            distanceToParticle = tol;
        }
    //int thisTmp=0;
    //int Bperiodic = boundaryVector[minIndex].periodic;
    //int BpointLine = boundaryVector[minIndex].pointLine;
    //int BsurfaceNumber = boundaryVector[minIndex].surfaceNumber;
    //int Bsurface = boundaryVector[minIndex].surface;
    //float Bx1 = boundaryVector[minIndex].x1;
    //float By1 = boundaryVector[minIndex].y1;
    //float Bz1 = boundaryVector[minIndex].z1;
    //float Bx2 = boundaryVector[minIndex].x2;
    //float By2 = boundaryVector[minIndex].y2;
    //float Bz2 = boundaryVector[minIndex].z2;
    //float Ba = boundaryVector[minIndex].a;
    //float Bb = boundaryVector[minIndex].b;
    //float Bc = boundaryVector[minIndex].c;
    //float Bd = boundaryVector[minIndex].d;
    //float Bplane_norm = boundaryVector[minIndex].plane_norm;
    //#if USE3DTETGEOM > 0
    //  float Bx3 = boundaryVector[minIndex].x3;
    //  float By3 = boundaryVector[minIndex].y3;
    //  float Bz3 = boundaryVector[minIndex].z3;
    //  float Barea = boundaryVector[minIndex].area;
    //#else
    //  float Bslope_dzdx = boundaryVector[minIndex].slope_dzdx;
    //  float Bintercept_z = boundaryVector[minIndex].intercept_z;
    //#endif     
    //float BZ = boundaryVector[minIndex].Z;
    //float Bamu = boundaryVector[minIndex].amu;
    //float Bpotential = boundaryVector[minIndex].potential;
	//
    //float BhitWall = boundaryVector[minIndex].hitWall;
    //float Blength = boundaryVector[minIndex].length;
    //float BdistanceToParticle = boundaryVector[minIndex].distanceToParticle;
    //float Bangle = boundaryVector[minIndex].angle;
    //float Bfd = boundaryVector[minIndex].fd;
    //float Bdensity = boundaryVector[minIndex].density;
    //float Bti = boundaryVector[minIndex].ti;
    //float Bne = boundaryVector[minIndex].ne;
    //float Bte = boundaryVector[minIndex].te;
    //float BdebyeLength = boundaryVector[minIndex].debyeLength;
    //float BlarmorRadius = boundaryVector[minIndex].larmorRadius;
    //   if(x0==0.0 && z > 1.0e-3 && minDistance<1.0e-9)
    //       thisTmp=1;
    }
    if (direction_type == 1)
    {
        if (boundaryVector[minIndex].slope_dzdx == 0)
        {
            directionUnitVector[0] = 0.0f;
            directionUnitVector[1] = 0.0f;
            directionUnitVector[2] = 1.0f * sgn(boundaryVector[minIndex].z1 - z);
        }
        else if (fabsf(boundaryVector[minIndex].slope_dzdx)>= 0.75f*tol)
        {
            
            directionUnitVector[0] = boundaryVector[minIndex].x1 - x;
            directionUnitVector[1] = 0.0f;
            directionUnitVector[2] = 0.0f;
        }
        else
        {
            directionUnitVector[0] = 1.0f * sgn((z - boundaryVector[minIndex].intercept_z)/(boundaryVector[minIndex].slope_dzdx) - x0);
            directionUnitVector[1] = 0.0f;
            directionUnitVector[2] = 1.0f * sgn(perp_dist)/(boundaryVector[minIndex].slope_dzdx);
        //std::cout << "sign boundarVec.slope  sign perp_dist " << sgn(boundaryVector[minIndex].slope_dzdx) << " " << sgn(perp_dist) << std::endl;
        }
        //std::cout << "direction_type 1 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0f;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
        //std::cout << "direction_type 2 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0f;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
        //std::cout << "direction_type 3 " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    }

    vectorMagnitude = sqrtf(directionUnitVector[0]*directionUnitVector[0] + directionUnitVector[1]*directionUnitVector[1]
                                + directionUnitVector[2]*directionUnitVector[2]);
    directionUnitVector[0] = directionUnitVector[0]/vectorMagnitude;
    directionUnitVector[1] = directionUnitVector[1]/vectorMagnitude;
    directionUnitVector[2] = directionUnitVector[2]/vectorMagnitude;
#endif   
#if BIASED_SURFACE > 0
    pot = boundaryVector[minIndex].potential;
    Emag = pot/(2.0f*boundaryVector[minIndex].ChildLangmuirDist)*expf(-minDistance/(2.0f*boundaryVector[minIndex].ChildLangmuirDist));
#else 
    angle = boundaryVector[minIndex].angle;    
    fd  =  0.98992f + 5.1220E-03f * angle  -
           7.0040E-04f  * powf(angle,2.0f) +
           3.3591E-05f  * powf(angle,3.0f) -
           8.2917E-07f  * powf(angle,4.0f) +
           9.5856E-09f   * powf(angle,5.0f) -
           4.2682E-11f  * powf(angle,6.0f);
    pot = boundaryVector[minIndex].potential;
     //std::cout << "potential and debye length " << pot << " " << boundaryVector[minIndex].debyeLength << " " << pot/boundaryVector[minIndex].debyeLength << std::endl;
    //std::cout << " larmorRad " << boundaryVector[minIndex].larmorRadius << std::endl;
        float debyeLength = boundaryVector[minIndex].debyeLength;
        float larmorRadius = boundaryVector[minIndex].larmorRadius;
        Emag = pot*(fd/(2.0f * boundaryVector[minIndex].debyeLength)*expf(-minDistance/(2.0f * boundaryVector[minIndex].debyeLength))+ (1.0f - fd)/(boundaryVector[minIndex].larmorRadius)*expf(-minDistance/boundaryVector[minIndex].larmorRadius) );
        float part1 = pot*(fd/(2.0f * boundaryVector[minIndex].debyeLength)*expf(-minDistance/(2.0f * boundaryVector[minIndex].debyeLength)));
        float part2 = pot*(1.0f - fd)/(boundaryVector[minIndex].larmorRadius)*expf(-minDistance/boundaryVector[minIndex].larmorRadius);
        //std::cout << "Emag " << Emag << std::endl;
        //std::cout << "fd " << fd << std::endl;
        //std::cout << "minDistance " << minDistance << std::endl;
#endif
    if(minDistance == 0.0f || boundaryVector[minIndex].larmorRadius == 0.0f)
    {
        Emag = 0.0f;
        directionUnitVector[0] = 0.0f;
        directionUnitVector[1] = 0.0f;
        directionUnitVector[2] = 0.0f;

    }
        Er = Emag*directionUnitVector[0];
        Et = Emag*directionUnitVector[1];
        E[2] = Emag*directionUnitVector[2];
        //std::cout << "Emag " << Emag << std::endl;
        //std::cout << "Min dist " << minDistance << std::endl;
        //std::cout << "r " << x << "z " << z << std::endl;
        //std::cout << "E components " << Er << " " << Et << " " << E[2] << std::endl;
        //std::cout << "direction unit vector " << directionUnitVector[0] << " " << directionUnitVector[1] << " " << directionUnitVector[2] << std::endl;
    
    //std::cout << "pos " << x << " " << y << " "<< z << " min Dist" << minDistance << "Efield " << Emag << std::endl;
#if USE3DTETGEOM > 0
            E[0] = Er;
            E[1] = Et;
#else
#if USECYLSYMM > 0
            //if cylindrical geometry
            float theta = atan2f(y,x0);
  
            E[0] = cosf(theta)*Er - sinf(theta)*Et;
            E[1] = sinf(theta)*Er + cosf(theta)*Et;
#else
            E[0] = Er;
            E[1] = Et;
#endif
#endif
            //std::cout << "Ex and Ey and Ez " << E[0] << " " << E[1] << " " << E[2] << std::endl;
   
      return minDistance;
}

struct move_boris { 
    Particles *particlesPointer;
    //int& tt;
    Boundary *boundaryVector;
    int nR_Bfield;
    int nZ_Bfield;
    float * BfieldGridRDevicePointer;
    float * BfieldGridZDevicePointer;
    float * BfieldRDevicePointer;
    float * BfieldZDevicePointer;
    float * BfieldTDevicePointer;
    int nR_Efield;
    int nY_Efield;
    int nZ_Efield;
    float * EfieldGridRDevicePointer;
    float * EfieldGridYDevicePointer;
    float * EfieldGridZDevicePointer;
    float * EfieldRDevicePointer;
    float * EfieldZDevicePointer;
    float * EfieldTDevicePointer;
    int nR_closeGeom_sheath;
    int nY_closeGeom_sheath;
    int nZ_closeGeom_sheath;
    int n_closeGeomElements_sheath;
    float* closeGeomGridr_sheath;
    float* closeGeomGridy_sheath;
    float* closeGeomGridz_sheath;
    int* closeGeom_sheath; 
    const float span;
    const int nLines;
    float magneticForce[3];
    float electricForce[3];
    move_boris(Particles *_particlesPointer, float _span, Boundary *_boundaryVector,int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,
            float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,
            float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer,
            int _nR_Efield,int _nY_Efield, int _nZ_Efield,
            float * _EfieldGridRDevicePointer,
            float * _EfieldGridYDevicePointer,
            float * _EfieldGridZDevicePointer,
            float * _EfieldRDevicePointer,
            float * _EfieldZDevicePointer,
            float * _EfieldTDevicePointer,
            int _nR_closeGeom, int _nY_closeGeom,int _nZ_closeGeom, int _n_closeGeomElements, float *_closeGeomGridr,float *_closeGeomGridy, float *_closeGeomGridz, int *_closeGeom)
        
        : particlesPointer(_particlesPointer), span(_span), boundaryVector(_boundaryVector), nLines(_nLines), nR_Bfield(_nR_Bfield), 
        nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), 
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), 
        BfieldTDevicePointer(_BfieldTDevicePointer),
        nR_Efield(_nR_Efield),nY_Efield(_nY_Efield), nZ_Efield(_nZ_Efield), 
        EfieldGridRDevicePointer(_EfieldGridRDevicePointer), 
        EfieldGridYDevicePointer(_EfieldGridYDevicePointer), 
        EfieldGridZDevicePointer(_EfieldGridZDevicePointer),
        EfieldRDevicePointer(_EfieldRDevicePointer), EfieldZDevicePointer(_EfieldZDevicePointer),
        EfieldTDevicePointer(_EfieldTDevicePointer),
   nR_closeGeom_sheath(_nR_closeGeom),nY_closeGeom_sheath(_nY_closeGeom), nZ_closeGeom_sheath(_nZ_closeGeom), n_closeGeomElements_sheath(_n_closeGeomElements), 
   closeGeomGridr_sheath(_closeGeomGridr),closeGeomGridy_sheath(_closeGeomGridy), closeGeomGridz_sheath(_closeGeomGridz), closeGeom_sheath(_closeGeom), magneticForce{0.0,0.0,0.0}, electricForce{0.0,0.0,0.0} {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx) { 
#ifdef __CUDACC__
#else
float initTime = 0.0f;
float interpETime = 0.0f;
float interpBTime = 0.0f;
float operationsTime = 0.0f;
#if USE_BOOST
cpu_timer timer;
cpu_times initTime0 = timer.elapsed();
#endif
#endif
            float v_minus[3]= {0.0f, 0.0f, 0.0f};
            float v_prime[3]= {0.0f, 0.0f, 0.0f};
            float position[3]= {0.0f, 0.0f, 0.0f};
	        float v[3]= {0.0f, 0.0f, 0.0f};
	        float E[3] = {0.0f, 0.0f, 0.0f};
#if USEPRESHEATHEFIELD > 0
            float PSE[3] = {0.0f, 0.0f, 0.0f};
#endif
	        float B[3] = {0.0f,0.0f,0.0f};
	        float dt = span;
	        float Bmag = 0.0f;
	        float q_prime = 0.0f;
            float coeff = 0.0f;
            int nSteps = floor( span / dt + 0.5f);
#if USESHEATHEFIELD > 0
            float minDist = 0.0f;
            int closestBoundaryIndex;
#endif
#if ODEINT ==	0 
        if(particlesPointer->hasLeaked[indx] == 0)
	{
	  if(particlesPointer->zprevious[indx] > particlesPointer->leakZ[indx])
	  {
	    particlesPointer->hasLeaked[indx] = 1;
	  }
	}
	        float qpE[3] = {0.0f,0.0f,0.0f};
	        float vmxB[3] = {0.0f,0.0f,0.0f};
	        float vpxB[3] = {0.0f,0.0f,0.0f};
	        float qp_vmxB[3] = {0.0f,0.0f,0.0f};
	        float c_vpxB[3] = {0.0f,0.0f,0.0f};
            vectorAssign(particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx],position);
            
            for ( int s=0; s<nSteps; s++ ) 
            {
#if USESHEATHEFIELD > 0
	          minDist = getE(particlesPointer->xprevious[indx], particlesPointer->yprevious[indx], particlesPointer->zprevious[indx],E,boundaryVector,nLines,nR_closeGeom_sheath,
                          nY_closeGeom_sheath,nZ_closeGeom_sheath,
                              n_closeGeomElements_sheath,closeGeomGridr_sheath,
                              closeGeomGridy_sheath,
                                   closeGeomGridz_sheath,closeGeom_sheath, closestBoundaryIndex);
              //std::cout << "Efield in boris " <<E[0] << " " << E[1] << " " <<  E[2] << std::endl;
              //std::cout << "Charge and Hitwall " << particlesPointer->charge[indx] << " " <<
               // particlesPointer->hitWall[indx]  << std::endl;
#endif

#if USEPRESHEATHEFIELD > 0
#if LC_INTERP==3
              
	        //float PSE2[3] = {0.0f, 0.0f, 0.0f};
                 interp3dVector(PSE,position[0], position[1], position[2],nR_Efield,nY_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridYDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer);
//E[0]= E[0] + PSE[0];
//E[1]= E[1] + PSE[1];
//E[2]= E[2] + PSE[2];
                 vectorAdd(E,PSE,E);
              //float a = interp3d(position[0], position[1], position[2],nR_Efield,nY_Efield,nZ_Efield,
                //                EfieldGridRDevicePointer,EfieldGridYDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer);
              //PSE[0] = 1.23;

#else
                 interp2dVector(&PSE[0],position[0], position[1], position[2],nR_Efield,nZ_Efield,
                     EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
                     EfieldZDevicePointer,EfieldTDevicePointer);
                 
                 vectorAdd(E,PSE,E);
#endif
#endif              
                interp2dVector(&B[0],position[0], position[1], position[2],nR_Bfield,nZ_Bfield,
                    BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                    BfieldZDevicePointer,BfieldTDevicePointer);        
                //std::cout << "Bfield and mass " <<B[0] << " " <<  B[1] <<" "<< B[2]<< " " << particlesPointer->amu[indx] << std::endl;    
                Bmag = vectorNorm(B);
	            q_prime = particlesPointer->charge[indx]*1.60217662e-19f/(particlesPointer->amu[indx]*1.6737236e-27f)*dt*0.5f;
                //std::cout << "charge, amu , dt " << particlesPointer->charge[indx] << " " << particlesPointer->amu[indx]<< " " << dt << std::endl;
                coeff = 2.0f*q_prime/(1.0f+(q_prime*Bmag)*(q_prime*Bmag));
                //std::cout << " Bmag " << Bmag << std::endl;
                //std::cout << " qprime coeff " << q_prime << " " << coeff << std::endl;
                vectorAssign(particlesPointer->vx[indx], particlesPointer->vy[indx], particlesPointer->vz[indx],v);
                //std::cout << "velocity " << v[0] << " " << v[1] << " " << v[2] << std::endl;
                //v_minus = v + q_prime*E;
               vectorScalarMult(q_prime,E,qpE);
               vectorAdd(v,qpE,v_minus);
               this->electricForce[0] = 2.0*qpE[0];
	       //std::cout << "e force " << q_prime << " " << PSE[0] << " " << PSE[1] << " " << PSE[2] << std::endl;
               this->electricForce[1] = 2.0*qpE[1];
               this->electricForce[2] = 2.0*qpE[2];
               //v_prime = v_minus + q_prime*(v_minus x B)
                vectorCrossProduct(v_minus,B,vmxB);
                vectorScalarMult(q_prime,vmxB,qp_vmxB);
                vectorAdd(v_minus,qp_vmxB,v_prime);       
               this->magneticForce[0] = qp_vmxB[0];
               this->magneticForce[1] = qp_vmxB[1];
               this->magneticForce[2] = qp_vmxB[2];
                
                //v = v_minus + coeff*(v_prime x B)
                vectorCrossProduct(v_prime, B, vpxB);
                vectorScalarMult(coeff,vpxB,c_vpxB);
                vectorAdd(v_minus, c_vpxB, v);
                
                //v = v + q_prime*E
                vectorAdd(v,qpE,v);
       //particlesPointer->test[indx] = Bmag; 
       //particlesPointer->test0[indx] = v[0]; 
       //particlesPointer->test1[indx] = v[1]; 
       //particlesPointer->test2[indx] = v[2];
        float ti_eV = 50.0;
	//std::cout << "ti dens tau_s " << ti_eV << " " << density << " " << tau_s << endl;
	float vTherm = sqrt(2*ti_eV*1.602e-19/particlesPointer->amu[indx]/1.66e-27);

      if(abs(v[2]) > vTherm)
      {
          v[2] = sgn(v[2])*vTherm;
          v[0] = 0.0;
          v[1] = 0.0;
      }
                float vxy00 = sqrt(vTherm*vTherm - v[2]*v[2]);
                float vxy01 = sqrt(v[1]*v[1]+ v[0]*v[0]);
		//std::cout << "vzNew vxy0 vxy " << vzNew << " " << vxy0 << " " << vxy << endl;
               v[0] = v[0]/vxy01*vxy00;///velocityCollisionsNorm; 
		       v[1] = v[1]/vxy01*vxy00;///velocityCollisionsNorm;


	       
           if(particlesPointer->hitWall[indx] == 0.0)
            {
                //std::cout << "updating r and v " << std::endl;
                particlesPointer->x[indx] = position[0] + v[0] * dt;
                particlesPointer->y[indx] = position[1] + v[1] * dt;
                particlesPointer->z[indx] = position[2] + v[2] * dt;
                particlesPointer->vx[indx] = v[0];
                particlesPointer->vy[indx] = v[1];
                particlesPointer->vz[indx] = v[2];    
              
//std::cout << "velocity " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    	    }
#endif

#if ODEINT == 1
        float m = p.amu*1.6737236e-27;
        float q_m = p.charge*1.60217662e-19/m;
        float r[3]= {0.0, 0.0, 0.0};
        float r2[3]= {0.0, 0.0, 0.0};
        float r3[3]= {0.0, 0.0, 0.0};
        float r4[3]= {0.0, 0.0, 0.0};
        float v2[3]= {0.0, 0.0, 0.0};
        float v3[3]= {0.0, 0.0, 0.0};
        float v4[3]= {0.0, 0.0, 0.0};
        float k1r[3]= {0.0, 0.0, 0.0};
        float k2r[3]= {0.0, 0.0, 0.0};
        float k3r[3]= {0.0, 0.0, 0.0};
        float k4r[3]= {0.0, 0.0, 0.0};
        float k1v[3]= {0.0, 0.0, 0.0};
        float k2v[3]= {0.0, 0.0, 0.0};
        float k3v[3]= {0.0, 0.0, 0.0};
        float k4v[3]= {0.0, 0.0, 0.0};
        float dtqm = dt*q_m;
        float vxB[3] = {0.0,0.0,0.0};
        float EplusvxB[3] = {0.0,0.0,0.0};
        float halfKr[3] = {0.0,0.0,0.0};
        float halfKv[3] = {0.0,0.0,0.0};
        float half = 0.5;
                v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;

                r[0] = p.xprevious;
                r[1] = p.yprevious;
	            r[2] = p.zprevious;
#ifdef __CUDACC__
#else
#if USE_BOOST
cpu_times initTime1 = timer.elapsed();
initTime = initTime + (initTime1.wall - initTime0.wall);
#endif
#endif
for ( int s=0; s<nSteps; s++ ) 
    {
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times operationsTime0 = timer.elapsed();
    cpu_times interpETime0 = timer.elapsed();
#endif
#endif
#if USESHEATHEFIELD > 0
    minDist = getE(r[0],r[1],r[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
          EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
          EfieldZDevicePointer,EfieldTDevicePointer);
                 
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times interpBTime0 = timer.elapsed();
#endif     
#endif
    interp2dVector(&B[0],r[0],r[1],r[2],nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
               BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k1r = dt*v
    vectorScalarMult(dt,v,k1r);
    /*
    k1r[0] = v[0]*dt;
    k1r[1] = v[1]*dt;
    k1r[2] = v[2]*dt;
    */
    //k1v = dt*q_m * (E + (v x B))
    vectorCrossProduct(v,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k1v);
    /*
    k1v[0] = dt*q_m*(E[0] + (v[1]*B[2] - v[2]*B[1]));
    k1v[1] = dt*q_m*(E[1] + (v[2]*B[0] - v[0]*B[2]));
    k1v[2] = dt*q_m*(E[2] + (v[0]*B[1] - v[1]*B[0]));
    */
    //r2 = r + 0.5*k1r
    vectorScalarMult(half,k1r,halfKr);
    vectorAdd(r,k1r,r2);
    /*
    r2[0] = r[0] + k1r[0]*0.5;
    r2[1] = r[1] + k1r[1]*0.5;
    r2[2] = r[2] + k1r[2]*0.5;
    */

    //v2 = v + 0.5*k1v
    vectorScalarMult(half,k1v,halfKv);
    vectorAdd(v, halfKv,v2);
        /*
    v2[0] = v[0] + k1v[0]*0.5;
    v2[1] = v[1] + k1v[1]*0.5;
    v2[2] = v[2] + k1v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0	  
    minDist = getE(r2[0],r2[1],r2[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif


    interp2dVector(&B[0],r2[0],r2[1],r2[2],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
             BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k2r = dt*v2
    vectorScalarMult(dt,v2,k2r);
    /*
    k2r[0] = v2[0]*dt;
    k2r[1] = v2[1]*dt;
    k2r[2] = v2[2]*dt;
    */
    //k2v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v2,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k2v);
    /*
    k2v[0] = dt*q_m*(E[0] + (v2[1]*B[2] - v2[2]*B[1]));
    k2v[1] = dt*q_m*(E[1] + (v2[2]*B[0] - v2[0]*B[2]));
    k2v[2] = dt*q_m*(E[2] + (v2[0]*B[1] - v2[1]*B[0]));
    */
    //r3 = r + 0.5*k2r
    vectorScalarMult(half,k2r,halfKr);
    vectorAdd(r,k2r,r3);
    /*
    r3[0] = r[0] + k2r[0]*0.5;
    r3[1] = r[1] + k2r[1]*0.5;
    r3[2] = r[2] + k2r[2]*0.5;
    */
    //v3 = v + 0.5*k2v
    vectorScalarMult(half,k2v,halfKv);
    vectorAdd(v, halfKv,v3);
    /*
    v3[0] = v[0] + k2v[0]*0.5;
    v3[1] = v[1] + k2v[1]*0.5;
    v3[2] = v[2] + k2v[2]*0.5;
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0	  
    minDist = getE(r3[0],r3[1],r3[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
    interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              

#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif
    interp2dVector(&B[0],r3[0],r3[1],r3[2],nR_Bfield,nZ_Bfield,
                 BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,
                 BfieldZDevicePointer,BfieldTDevicePointer);        
                
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif
    //k3r = dt*v3
    vectorScalarMult(dt,v3,k3r);
    /*
    k3r[0] = v3[0]*dt;
    k3r[1] = v3[1]*dt;
    k3r[2] = v3[2]*dt;
    */
    //k3v = dt*qm*(E + (v x B))
    vectorCrossProduct(v3,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k3v);
    /*
    k3v[0] = dt*q_m*(E[0] + (v3[1]*B[2] - v3[2]*B[1]));
    k3v[1] = dt*q_m*(E[1] + (v3[2]*B[0] - v3[0]*B[2]));
    k3v[2] = dt*q_m*(E[2] + (v3[0]*B[1] - v3[1]*B[0]));
    */
    //r4 = r + k3r
    vectorAdd(r, k3r,r4);
    /*
    r4[0] = r[0] + k3r[0];
    r4[1] = r[1] + k3r[1];
    r4[2] = r[2] + k3r[2];
    */
    //v4 = v + k3v
    vectorAdd(v, k3v, v4);
        /*
    v4[0] = v[0] + k3v[0];
    v4[1] = v[1] + k3v[1];
    v4[2] = v[2] + k3v[2];
    */
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpETime0 = timer.elapsed();
#endif
#endif

#if USESHEATHEFIELD > 0            
	minDist = getE(r4[0],r4[1],r4[2],E,boundaryVector,nLines);
#endif
#if USEPRESHEATHEFIELD > 0
   interp2dVector(&PSE[0],p.xprevious,p.yprevious,p.zprevious,nR_Efield,nZ_Efield,
               EfieldGridRDevicePointer,EfieldGridZDevicePointer,EfieldRDevicePointer,
               EfieldZDevicePointer,EfieldTDevicePointer);
    vectorAdd(E,PSE,E);
#endif              
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime0 = timer.elapsed();
#endif
#endif

    interp2dVector(&B[0],r4[0],r4[1],r4[2],nR_Bfield,nZ_Bfield,
                        BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                        BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);        
#ifdef __CUDACC__
#else
#if USE_BOOST
    interpBTime1 = timer.elapsed();
    interpETime = interpETime + (interpBTime0.wall - interpETime0.wall);
    interpBTime = interpBTime + (interpBTime1.wall - interpBTime0.wall);
#endif
#endif

    //k4r = dt*v4
    vectorScalarMult(dt,v4,k4r);
    /*
   k4r[0] = v4[0]*dt;
   k4r[1] = v4[1]*dt;
   k4r[2] = v4[2]*dt;
   */
    //k4v = dt*q_m*(E + (v x B))
    vectorCrossProduct(v4,B,vxB);
    vectorAdd(E,vxB,EplusvxB);
    vectorScalarMult(dtqm,EplusvxB,k4v);
    /*
   k4v[0] = dt*q_m*(E[0] + (v4[1]*B[2] - v4[2]*B[1]));
   k4v[1] = dt*q_m*(E[1] + (v4[2]*B[0] - v4[0]*B[2]));
   k4v[2] = dt*q_m*(E[2] + (v4[0]*B[1] - v4[1]*B[0]));
   */
   p.x = r[0] + (k1r[0] + 2*k2r[0] + 2*k3r[0] + k4r[0])/6;
   p.y = r[1] + (k1r[1] + 2*k2r[1] + 2*k3r[1] + k4r[1])/6;
   p.z = r[2] + (k1r[2] + 2*k2r[2] + 2*k3r[2] + k4r[2])/6;
   p.vx = v[0] + (k1v[0] + 2*k2v[0] + 2*k3v[0] + k4v[0])/6;
   p.vy = v[1] + (k1v[1] + 2*k2v[1] + 2*k3v[1] + k4v[1])/6;
   p.vz = v[2] + (k1v[2] + 2*k2v[2] + 2*k3v[2] + k4v[2])/6;
#ifdef __CUDACC__
#else
#if USE_BOOST
    cpu_times operationsTime1 = timer.elapsed();
    operationsTime = operationsTime1.wall - operationsTime0.wall;
#endif
#endif
//std::cout << "Operations Time: " << operationsTime <<std::endl;
//std::cout << "Efield Interpolation Time: " << interpETime <<std::endl;
//std::cout << "Bfield Interpolation Time: " << interpBTime <<std::endl;
//std::cout << "Init Time: " << initTime <<std::endl;
            }
#endif
        }
    } 
};

#endif
