#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include "BoundaryModifiable.h"
#include <math.h>

/*template <typename T>
CUDA_CALLABLE_MEMBER_DEVICE
int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}*/

struct geometry_check { 
    Particles *particlesPointer;
    BoundaryModifiable *boundaryMod;
    const int nLines;
    Boundary  *boundaryVector;
    float dt;
    int tt;
    int nR_closeGeom;
    int nY_closeGeom;
    int nZ_closeGeom;
    int n_closeGeomElements;
    float* closeGeomGridr;
    float* closeGeomGridy;
    float* closeGeomGridz;
    int* closeGeom;

    geometry_check(Particles *_particlesPointer, BoundaryModifiable *_boundaryMod, int _nLines,Boundary * _boundaryVector, float _dt, int _tt, int _nR_closeGeom, int _nY_closeGeom, int _nZ_closeGeom, int _n_closeGeomElements, float *_closeGeomGridr, float *_closeGeomGridy, float *_closeGeomGridz, int *_closeGeom) : 
        particlesPointer(_particlesPointer), boundaryMod(_boundaryMod), nLines(_nLines), boundaryVector(_boundaryVector), dt(_dt), tt(_tt), nR_closeGeom(_nR_closeGeom), nY_closeGeom(_nY_closeGeom), nZ_closeGeom(_nZ_closeGeom), n_closeGeomElements(_n_closeGeomElements), closeGeomGridr(_closeGeomGridr), closeGeomGridy(_closeGeomGridy), closeGeomGridz(_closeGeomGridz), closeGeom(_closeGeom) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 
    //std::cout << "geometry check particle x" << particlesPointer->x[indx] << particlesPointer->x[indx]previous <<std::endl;
    //std::cout << "geometry check particle y" << particlesPointer->y[indx] << particlesPointer->y[indx]previous <<std::endl;
    //std::cout << "geometry check particle z" << particlesPointer->z[indx] << particlesPointer->z[indx]previous <<std::endl;
    //std::cout << "geometry check particle hitwall" << p.hitWall <<std::endl;
	    if(particlesPointer->hitWall[indx] == 0.0)
        { 
#if USE3DTETGEOM > 0
            float p0[3] = {particlesPointer->xprevious[indx],
                           particlesPointer->yprevious[indx],
                           particlesPointer->zprevious[indx]};
            float p1[3] = {particlesPointer->x[indx],
                           particlesPointer->y[indx],
                           particlesPointer->z[indx]};

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
            float normalVector[3] = {0.0,0.0,0.0};
            float crossABAp[3] = {0.0,0.0,0.0};
            float crossBCBp[3] = {0.0,0.0,0.0};
            float crossCACp[3] = {0.0,0.0,0.0};
            float signDot0 = 0.0;
            float signDot1 = 0.0;
            float signDot2 = 0.0;
            float totalSigns = 0.0;
            int nBoundariesCrossed = 0;
            int boundariesCrossed[6] = {0,0,0,0,0,0};
            /*
            if(particlesPointer->xprevious[indx] < 0 || particlesPointer->yprevious[indx] < 0 ||
                    particlesPointer->zprevious[indx]< 0)
            {
            std::cout << "pos " << particlesPointer->xprevious[indx] << " " 
                << particlesPointer->yprevious[indx]
                << " " << particlesPointer->zprevious[indx]  << std::endl;
            }
             */
#if GEOM_HASH == 1
            float r_position = particlesPointer->xprevious[indx];
            float dr = closeGeomGridr[1] - closeGeomGridr[0];
            float dz = closeGeomGridz[1] - closeGeomGridz[0];    
            int rInd = floor((r_position - closeGeomGridr[0])/dr + 0.5f);
            int zInd = floor((particlesPointer->zprevious[indx] - closeGeomGridz[0])/dz + 0.5f);
            int i=0;
            float dy = closeGeomGridy[1] - closeGeomGridy[0];    
            int yInd = floor((particlesPointer->yprevious[indx] - closeGeomGridy[0])/dy + 0.5f);
            if(rInd < 0 || yInd < 0 || zInd < 0)
            {
            //std::cout << "indices " << rInd << " " << yInd << " " << zInd << std::endl;
            //std::cout << "pos " << r_position << " " << particlesPointer->yprevious[indx]
            //    << " " << particlesPointer->zprevious[indx]  << std::endl;
              i = closeGeom[-1];
            }
            for (int j=0; j< n_closeGeomElements; j++)
            {
                i = closeGeom[zInd*nY_closeGeom*nR_closeGeom*n_closeGeomElements
                            + yInd*nR_closeGeom*n_closeGeomElements + rInd*n_closeGeomElements + j]; 
#else
            for (int i=0; i<nLines; i++)
                    {
#endif
                      a = boundaryVector[i].a;
                      b = boundaryVector[i].b;
                      c = boundaryVector[i].c;
                      d = boundaryVector[i].d;
                      plane_norm = boundaryVector[i].plane_norm;
                      pointToPlaneDistance0 = (a*p0[0] + b*p0[1] + c*p0[2] + d)/plane_norm;
                      pointToPlaneDistance1 = (a*p1[0] + b*p1[1] + c*p1[2] + d)/plane_norm;
                    
                      signPoint0 = sgn(pointToPlaneDistance0);
                      signPoint1 = sgn(pointToPlaneDistance1);

                      if (signPoint0 != signPoint1)
                      {
                        t = -(a*p0[0] + b*p0[1] + c*p0[2] + d)/
                            (a*(p1[0] - p0[0]) + b*(p1[1] - p0[1]) + c*(p1[2] - p0[2]));
                        vectorAssign(p0[0] + t*(p1[0] - p0[0]), 
                                p0[1] + t*(p1[1] - p0[1]),
                                p0[2] + t*(p1[2] - p0[2]), p);

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

                        signDot0 = sgn(vectorDotProduct(crossABAp,normalVector));
                        signDot1 = sgn(vectorDotProduct(crossBCBp,normalVector));
                        signDot2 = sgn(vectorDotProduct(crossCACp,normalVector));
                        totalSigns = abs(signDot0 + signDot1 + signDot2);

                        if (totalSigns == 3.0)
                        {
                            boundariesCrossed[nBoundariesCrossed] = i; 
                            nBoundariesCrossed++;
                            particlesPointer->hitWall[indx] = 1.0;
                            particlesPointer->xprevious[indx] = p[0];
                            particlesPointer->yprevious[indx] = p[1];
                            particlesPointer->zprevious[indx] = p[2];
                            particlesPointer->x[indx] = p[0];
                            particlesPointer->y[indx] = p[1];
                            particlesPointer->z[indx] = p[2];
                            particlesPointer->wallHit[indx] = i; 
#if USE_CUDA > 0
                            atomicAdd(&boundaryVector[i].impacts, 1.0);
#else
                            //boundaryMod->impacts[i] = boundaryMod->impacts[i] + 1.0;
                            boundaryVector[i].impacts = boundaryVector[i].impacts +  1.0;
#endif
                        }   
                     }
                        if (nBoundariesCrossed == 0)
                        {    
                            particlesPointer->xprevious[indx] = particlesPointer->x[indx];
                            particlesPointer->yprevious[indx] = particlesPointer->y[indx];
                            particlesPointer->zprevious[indx] = particlesPointer->z[indx];
                        }
                    } 
#else            
#if USECYLSYMM > 0
           float pdim1 = sqrt(particlesPointer->x[indx]*particlesPointer->x[indx] + particlesPointer->y[indx]*particlesPointer->y[indx]);
           float pdim1previous = sqrt(particlesPointer->xprevious[indx]*particlesPointer->xprevious[indx] + particlesPointer->yprevious[indx]*particlesPointer->yprevious[indx]); 
           float theta0 = atanf(particlesPointer->yprevious[indx]/particlesPointer->xprevious[indx]);
           float theta1 = atanf(particlesPointer->y[indx]/particlesPointer->x[indx]);
           float thetaNew = 0;
#else      
           float pdim1 = particlesPointer->x[indx];
           float pdim1previous = particlesPointer->xprevious[indx];
#endif            
            float particle_slope = (particlesPointer->z[indx] - particlesPointer->zprevious[indx])/(pdim1 - pdim1previous);
            float particle_intercept = -particle_slope*pdim1 + particlesPointer->z[indx];
            float intersectionx[2] = {};
            //intersectionx = new float[nPoints];
            float intersectiony[2] = {};
            //intersectiony = new float[nPoints];
            float distances[2] = {};
            //distances = new float[nPoints];
            int intersectionIndices[2] = {};
            float tol_small = 1e-12f;       
            float tol = 1e12f;
	        int nIntersections = 0;
            float signPoint;
            float signPoint0;
            float signLine1;
            float signLine2;
            float minDist = 1e12f;
            int minDistInd = 0;
                //std::cout << "particle slope " << particle_slope << " " << particle_intercept << std::endl;
                            //std::cout << "r " << boundaryVector[0].x1 << " " << boundaryVector[0].x1 << " " << boundaryVector[0].slope_dzdx << std::endl;
                            //std::cout << "r0 " << particlesPointer->x[indx]previous << " " << particlesPointer->y[indx]previous << " " << particlesPointer->z[indx]previous<< std::endl;
#if GEOM_HASH == 1
#if USECYLSYMM > 0
            float r_position = sqrtf(particlesPointer->xprevious[indx]*particlesPointer->xprevious[indx] + particlesPointer->yprevious[indx]*particlesPointer->yprevious[indx]);
#else   
            float r_position = particlesPointer->xprevious[indx];
#endif
            float dr = closeGeomGridr[1] - closeGeomGridr[0];
            float dz = closeGeomGridz[1] - closeGeomGridz[0];    
            int rInd = floor((r_position - closeGeomGridr[0])/dr + 0.5f);
            int zInd = floor((particlesPointer->zprevious[indx] - closeGeomGridz[0])/dz + 0.5f);
            int i;
            for (int j=0; j< n_closeGeomElements; j++)
            {
                i = closeGeom[zInd*nR_closeGeom*n_closeGeomElements + rInd*n_closeGeomElements + j];

#else
            for (int i=0; i<nLines; i++)
            {
#endif    
               //std::cout << "vert geom " << i << "  " << fabs(boundaryVector[i].slope_dzdx) << " " << tol << std::endl;
                    if (fabsf(boundaryVector[i].slope_dzdx) >= tol*0.75f)
                    {
                        signPoint = sgn(pdim1 - boundaryVector[i].x1);
                        signPoint0 = sgn(pdim1previous - boundaryVector[i].x1);
                        //std::cout << "signpoint1 " << signPoint << " " << signPoint0 << std::endl;
                    }
                    else
                    {
                        signPoint = sgn(particlesPointer->z[indx] - pdim1*boundaryVector[i].slope_dzdx - boundaryVector[i].intercept_z);
                        signPoint0 = sgn(particlesPointer->zprevious[indx] - pdim1previous*boundaryVector[i].slope_dzdx - boundaryVector[i].intercept_z);
                        //std::cout << "signpoint2 " << signPoint << " " << signPoint0 << std::endl;
                    }

                    if (signPoint != signPoint0)
                    {
                       if (fabsf(particle_slope)>= tol*0.75f)
                               {
                                  // std::cout << " isinf catch " << std::endl;
                                particle_slope = tol;
                               } 
                       if (fabsf(particle_slope) >= tol*0.75f)
                       {
                           signLine1 = sgn(boundaryVector[i].x1 - pdim1);
                           signLine2 = sgn(boundaryVector[i].x2 - pdim1);
                          // std::cout << "signlines3 " << signLine1 << " " << signLine2 << std::endl;
                       }
                       else
                       {
                        signLine1 = sgn(boundaryVector[i].z1 - boundaryVector[i].x1*particle_slope - particle_intercept);
                        signLine2 = sgn(boundaryVector[i].z2 - boundaryVector[i].x2*particle_slope - particle_intercept);
                       }
                        //std::cout << "signLines " << signLine1 << " " << signLine2 << std::endl;
                        //std::cout << "bound vec points " << boundaryVector[i].z1 << " " << boundaryVector[i].x1 << 
                          // " " << boundaryVector[i].z2 << " " << boundaryVector[i].x2 << std::endl; 
                        if (signLine1 != signLine2)
                        {
                            intersectionIndices[nIntersections] = i;
                            nIntersections++;

                            //std::cout << "nintersections " << nIntersections << std::endl;
                           // std::cout << fabs(particlesPointer->x[indx] - particlesPointer->x[indx]previous) << tol_small << std::endl;        
                            if (fabsf(pdim1 - pdim1previous) < tol_small)
                            {
                              //  std::cout << "vertical line" << std::cout;
                                intersectionx[nIntersections-1] = pdim1previous;
                                intersectiony[nIntersections-1] = intersectionx[nIntersections-1]*boundaryVector[i].slope_dzdx +
                                                                    boundaryVector[i].intercept_z;
                            }
                            else
                            {
                               //std::cout << "not vertical line" << std::endl;
                            //std::cout << 0.0*7.0 << " " << i << " " << nParam << " " << lines[i*nParam+4] << "  " <<tol << std::endl;
                               //std::cout << "boundaryVector slope " << boundaryVector[i].slope_dzdx << " " << tol*0.75 <<std::endl; 
                               if (fabsf(boundaryVector[i].slope_dzdx) >= tol*0.75f)
                                {
                                    intersectionx[nIntersections-1] = boundaryVector[i].x1;
                                }
                                else
                                {
                                    intersectionx[nIntersections-1] = (boundaryVector[i].intercept_z - particle_intercept)/
                                                (particle_slope - boundaryVector[i].slope_dzdx);
                                  //  std::cout << "in this else "<< intersectionx[nIntersections -1] << std::endl;
                                }
                                intersectiony[nIntersections-1] = intersectionx[nIntersections-1]*particle_slope
                                                                               + particle_intercept;
                            }
                        }
                    }
                }
       //if(particlesPointer->hitWall[indx] == 0.0)
         // {
            if (nIntersections ==0) 
            {
                particlesPointer->xprevious[indx] = particlesPointer->x[indx];
                particlesPointer->yprevious[indx] = particlesPointer->y[indx];
                particlesPointer->zprevious[indx] = particlesPointer->z[indx];

                //std::cout << "r " << particlesPointer->x[indx] << " " << particlesPointer->y[indx] << " " << particlesPointer->z[indx] << std::endl;
                //std::cout << "r0 " << particlesPointer->x[indx]previous << " " << particlesPointer->y[indx]previous << " " << particlesPointer->z[indx]previous<< std::endl;
            }
            else if (nIntersections ==1)
            {
                particlesPointer->hitWall[indx] = 1.0f;
                particlesPointer->wallIndex[indx] = intersectionIndices[0];
                if (particle_slope >= tol*0.75f)
                {
#if USECYLSYMM > 0
                    thetaNew = theta0 + (intersectiony[0] - particlesPointer->zprevious[indx])/(particlesPointer->z[indx] - particlesPointer->zprevious[indx])*(theta1 - theta0);
                   particlesPointer->y[indx] = intersectionx[0]*sinf(thetaNew);
#else                    
                    particlesPointer->y[indx] = particlesPointer->yprevious[indx] + (intersectiony[0] - particlesPointer->zprevious[indx])/(particlesPointer->z[indx] - particlesPointer->zprevious[indx])*(particlesPointer->y[indx] - particlesPointer->yprevious[indx]);
#endif     
                }
                else
                {
#if USECYLSYMM > 0
                thetaNew = theta0 + (intersectionx[0] - pdim1previous)/(pdim1 - pdim1previous)*(theta1 - theta0);    
                particlesPointer->y[indx] = intersectionx[0]*sinf(thetaNew);
#else                    
                particlesPointer->y[indx] = particlesPointer->yprevious[indx] + (intersectionx[0] - particlesPointer->xprevious[indx])/(particlesPointer->x[indx] - particlesPointer->xprevious[indx])*(particlesPointer->y[indx] - particlesPointer->yprevious[indx]);
#endif                
                }
#if USECYLSYMM > 0
                particlesPointer->x[indx] = intersectionx[0]*cosf(thetaNew);
#else                
                particlesPointer->x[indx] = intersectionx[0];
#endif     
                particlesPointer->z[indx] = intersectiony[0];
                //std::cout << "nInt = 1 position " << intersectionx[0] << " " << intersectiony[0]  << std::endl;
            }
            else
            {
                //std::cout << "nInts greater than 1 " << nIntersections << std::endl;
                for (int i=0; i<nIntersections; i++)
                {
                    distances[i] = (pdim1previous - intersectionx[i])*(pdim1previous - intersectionx[i]) + 
                        (particlesPointer->zprevious[indx] - intersectiony[i])*(particlesPointer->zprevious[indx] - intersectiony[i]);
                    if (distances[i] < minDist)
                    {
                        minDist = distances[i];
                        minDistInd = i;
                    }
                }

                particlesPointer->wallIndex[indx] = intersectionIndices[minDistInd];
                particlesPointer->hitWall[indx] = 1.0f;
#if USECYLSYMM > 0
       thetaNew = theta0 + (intersectionx[minDistInd] - pdim1previous)/(pdim1 - pdim1previous)*(theta1 - theta0);
       particlesPointer->y[indx] = intersectionx[minDistInd]*cosf(thetaNew);
       particlesPointer->x[indx] = intersectionx[minDistInd]*cosf(thetaNew);
#else
                particlesPointer->y[indx] = particlesPointer->yprevious[indx] + (intersectionx[minDistInd] - pdim1previous)/(pdim1 - pdim1previous)*(particlesPointer->y[indx] - particlesPointer->yprevious[indx]);
                particlesPointer->x[indx] = intersectionx[minDistInd];
#endif                
                particlesPointer->z[indx] = intersectiony[minDistInd];
            }

#if USECYLSYMM > 0 
#else            
            if (boundaryVector[nLines].periodic)
            {
                if (particlesPointer->y[indx] < boundaryVector[nLines].y1)
                {
                    particlesPointer->y[indx] = boundaryVector[nLines].y2;//  - (boundaryVector[nLines].y1 - particlesPointer->y[indx]);
                }
                else if (particlesPointer->y[indx] > boundaryVector[nLines].y2)
                {
                    particlesPointer->y[indx] = boundaryVector[nLines].y1;//  + (particlesPointer->y[indx] - boundaryVector[nLines].y2);
                }
            }
            else
            {
                if (particlesPointer->y[indx] < boundaryVector[nLines].y1)
                {
                    particlesPointer->hitWall[indx] = 1.0f;
                }
                else if (particlesPointer->y[indx] > boundaryVector[nLines].y2)
                {
                    particlesPointer->hitWall[indx] = 1.0f;
                }
            }
#endif
#endif        
            if (particlesPointer->hitWall[indx] == 1.0)
            {
                particlesPointer->transitTime[indx] = tt*dt;
            }    
        }

    //std::cout << "2geometry check particle x" << particlesPointer->x[indx] << particlesPointer->x[indx]previous <<std::endl;
    //std::cout << "2geometry check particle y" << particlesPointer->y[indx] << particlesPointer->y[indx]previous <<std::endl;
    //std::cout << "2geometry check particle z" << particlesPointer->z[indx] << particlesPointer->z[indx]previous <<std::endl;
    }
};

#endif
