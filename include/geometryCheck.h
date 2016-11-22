#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <math.h>

/*template <typename T>
CUDA_CALLABLE_MEMBER_DEVICE
int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}*/

struct geometry_check { 
    Particles *particlesPointer;
    const int nLines;
    Boundary * boundaryVector;
    float dt;
    int tt;

    geometry_check(Particles *_particlesPointer, int _nLines,Boundary * _boundaryVector, float _dt, int _tt) : 
        particlesPointer(_particlesPointer), nLines(_nLines), boundaryVector(_boundaryVector), dt(_dt), tt(_tt) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 
    //std::cout << "geometry check particle x" << particlesPointer->x[indx] << particlesPointer->x[indx]previous <<std::endl;
    //std::cout << "geometry check particle y" << particlesPointer->y[indx] << particlesPointer->y[indx]previous <<std::endl;
    //std::cout << "geometry check particle z" << particlesPointer->z[indx] << particlesPointer->z[indx]previous <<std::endl;
    //std::cout << "geometry check particle hitwall" << p.hitWall <<std::endl;
	    if(particlesPointer->hitWall[indx] == 0.0)
        {  
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
            for (int i=0; i<nLines; i++)
                {   //std::cout << "vert geom " << i << "  " << fabs(boundaryVector[i].slope_dzdx) << " " << tol << std::endl;
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
            if (particlesPointer->hitWall[indx] == 1)
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
