#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particle.h"
#include "Boundary.h"
#include <math.h>

/*template <typename T>
CUDA_CALLABLE_MEMBER_DEVICE
int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}*/

struct geometry_check { 

    const int nLines;
    Boundary * boundaryVector;

    geometry_check(int _nLines,Boundary * _boundaryVector) : nLines(_nLines), boundaryVector(_boundaryVector) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 
    //std::cout << "geometry check particle x" << p.x << p.xprevious <<std::endl;
    //std::cout << "geometry check particle y" << p.y << p.yprevious <<std::endl;
    //std::cout << "geometry check particle z" << p.z << p.zprevious <<std::endl;
    //std::cout << "geometry check particle hitwall" << p.hitWall <<std::endl;
	    if(p.hitWall == 0.0)
        {  
#if USECYLSYMM > 0
           float pdim1 = sqrt(p.x*p.x + p.y*p.y);
           float pdim1previous = sqrt(p.xprevious*p.xprevious + p.yprevious*p.yprevious); 
           float theta0 = atan(p.yprevious/p.xprevious);
           float theta1 = atan(p.y/p.x);
           float thetaNew = 0;
#else      
           float pdim1 = p.x;
           float pdim1previous = p.xprevious;
#endif            
            float particle_slope = (p.z - p.zprevious)/(pdim1 - pdim1previous);
            float particle_intercept = -particle_slope*pdim1 + p.z;
            float intersectionx[2] = {};
            //intersectionx = new float[nPoints];
            float intersectiony[2] = {};
            //intersectiony = new float[nPoints];
            float distances[2] = {};
            //distances = new float[nPoints];
            int intersectionIndices[2] = {};
            float tol_small = 1e-12;       
            float tol = 1e12;
	        int nIntersections = 0;
            float signPoint;
            float signPoint0;
            float signLine1;
            float signLine2;
            float minDist = 1e12;
            int minDistInd = 0;
                //std::cout << "particle slope " << particle_slope << " " << particle_intercept << std::endl;
                            //std::cout << "r " << boundaryVector[0].x1 << " " << boundaryVector[0].x1 << " " << boundaryVector[0].slope_dzdx << std::endl;
                            //std::cout << "r0 " << p.xprevious << " " << p.yprevious << " " << p.zprevious<< std::endl;
            for (int i=0; i<nLines; i++)
                {   //std::cout << "vert geom " << i << "  " << fabs(boundaryVector[i].slope_dzdx) << " " << tol << std::endl;
                    if (fabs(boundaryVector[i].slope_dzdx) >= tol*0.75)
                    {
                        signPoint = sgn(pdim1 - boundaryVector[i].x1);
                        signPoint0 = sgn(pdim1previous - boundaryVector[i].x1);
                        //std::cout << "signpoint1 " << signPoint << " " << signPoint0 << std::endl;
                    }
                    else
                    {
                        signPoint = sgn(p.z - pdim1*boundaryVector[i].slope_dzdx - boundaryVector[i].intercept_z);
                        signPoint0 = sgn(p.zprevious - pdim1previous*boundaryVector[i].slope_dzdx - boundaryVector[i].intercept_z);
                        //std::cout << "signpoint2 " << signPoint << " " << signPoint0 << std::endl;
                    }

                    if (signPoint != signPoint0)
                    {
                       if (fabs(particle_slope)>= tol*0.75)
                               {
                                  // std::cout << " isinf catch " << std::endl;
                                particle_slope = tol;
                               } 
                       if (fabs(particle_slope) >= tol*0.75)
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
                           // std::cout << fabs(p.x - p.xprevious) << tol_small << std::endl;        
                            if (fabs(pdim1 - pdim1previous) < tol_small)
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
                               if (fabs(boundaryVector[i].slope_dzdx) >= tol*0.75)
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
                p.xprevious = p.x;
                p.yprevious = p.y;
                p.zprevious = p.z;

                //std::cout << "r " << p.x << " " << p.y << " " << p.z << std::endl;
                //std::cout << "r0 " << p.xprevious << " " << p.yprevious << " " << p.zprevious<< std::endl;
            }
            else if (nIntersections ==1)
            {
                p.hitWall = 1.0;
                p.wallIndex = intersectionIndices[0];
                if (particle_slope >= tol*0.75)
                {
#if USECYLSYMM > 0
                    thetaNew = theta0 + (intersectiony[0] - p.zprevious)/(p.z - p.zprevious)*(theta1 - theta0);
                   p.y = intersectionx[0]*sin(thetaNew);
#else                    
                    p.y = p.yprevious + (intersectiony[0] - p.zprevious)/(p.z - p.zprevious)*(p.y - p.yprevious);
#endif     
                }
                else
                {
#if USECYLSYMM > 0
                thetaNew = theta0 + (intersectionx[0] - pdim1previous)/(pdim1 - pdim1previous)*(theta1 - theta0);    
                p.y = intersectionx[0]*sin(thetaNew);
#else                    
                p.y = p.yprevious + (intersectionx[0] - p.xprevious)/(p.x - p.xprevious)*(p.y - p.yprevious);
#endif                
                }
#if USECYLSYMM > 0
                p.x = intersectionx[0]*cos(thetaNew);
#else                
                p.x = intersectionx[0];
#endif     
                p.z = intersectiony[0];
                //std::cout << "nInt = 1 position " << intersectionx[0] << " " << intersectiony[0]  << std::endl;
            }
            else
            {
                //std::cout << "nInts greater than 1 " << nIntersections << std::endl;
                for (int i=0; i<nIntersections; i++)
                {
                    distances[i] = (pdim1previous - intersectionx[i])*(pdim1previous - intersectionx[i]) + 
                        (p.zprevious - intersectiony[i])*(p.zprevious - intersectiony[i]);
                    if (distances[i] < minDist)
                    {
                        minDist = distances[i];
                        minDistInd = i;
                    }
                }

                p.wallIndex = intersectionIndices[minDistInd];
                p.hitWall = 1.0;
#if USECYLSYMM > 0
       thetaNew = theta0 + (intersectionx[minDistInd] - pdim1previous)/(pdim1 - pdim1previous)*(theta1 - theta0);
       p.y = intersectionx[minDistInd]*cos(thetaNew);
       p.x = intersectionx[minDistInd]*cos(thetaNew);
#else
                p.y = p.yprevious + (intersectionx[minDistInd] - pdim1previous)/(pdim1 - pdim1previous)*(p.y - p.yprevious);
                p.x = intersectionx[minDistInd];
#endif                
                p.z = intersectiony[minDistInd];
            }

#if USECYLSYMM > 0 
#else            
            if (boundaryVector[nLines].periodic)
            {
                if (p.y < boundaryVector[nLines].y1)
                {
                    p.y = boundaryVector[nLines].y2;//  - (boundaryVector[nLines].y1 - p.y);
                }
                else if (p.y > boundaryVector[nLines].y2)
                {
                    p.y = boundaryVector[nLines].y1;//  + (p.y - boundaryVector[nLines].y2);
                }
            }
            else
            {
                if (p.y < boundaryVector[nLines].y1)
                {
                    p.hitWall = 1.0;
                }
                else if (p.y > boundaryVector[nLines].y2)
                {
                    p.hitWall = 1.0;
                }
            }
#endif        
        }

    //std::cout << "2geometry check particle x" << p.x << p.xprevious <<std::endl;
    //std::cout << "2geometry check particle y" << p.y << p.yprevious <<std::endl;
    //std::cout << "2geometry check particle z" << p.z << p.zprevious <<std::endl;
        }
};

#endif
