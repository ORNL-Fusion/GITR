#ifndef _GEOM_
#define _GEOM_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particle.h"

template <typename T>
CUDA_CALLABLE_MEMBER_DEVICE
int sgn(T val) {
        return (T(0) < val) - (val < T(0));
}

struct geometry_check { 

    const double var;

    geometry_check(double _var) : var(_var) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {   int nPoints = 4;
            double lines[28] =    { 0.0173, 0.0300, -0.0173, -0.0300, 1.7321, 0.0,0.0693, 
                                 -0.0173, -0.0300, -0.0200, -0.0300, 0.0, -0.0300, 0.0027,
                                 -0.0200, -0.0300, -0.0200, 0.0300, 1e12, 1e12,0.0600,
                                -0.0200,0.0300,0.0173,0.0300,0.0,0.0300,0.0373 };
            int nParam = 7;
            double particle_slope = (p.z - p.zprevious)/(p.x - p.xprevious);
            double particle_intercept = -particle_slope*p.x + p.z;
            double intersectionx[4] = {0.0, 0.0, 0.0, 0.0};
            //intersectionx = new double[nPoints];
            double intersectiony[4] = {0.0, 0.0, 0.0, 0.0};
            //intersectiony = new double[nPoints];
            double distances[4] = {0.0, 0.0, 0.0, 0.0};
            //distances = new double[nPoints];
            double tol_small = 1e-12;       
            double tol = 1e12;
	        int nIntersections = 0;
            double signPoint;
            double signPoint0;
            double signLine1;
            double signLine2;
            double minDist = 1e12;
            int minDistInd = 0;

                            //std::cout << "r " << p.x << " " << p.y << " " << p.z << std::endl;
                            //std::cout << "r0 " << p.xprevious << " " << p.yprevious << " " << p.zprevious<< std::endl;
            for (int i=0; i<nPoints; i++)
                {
                    if (fabs(lines[i*nParam+4]) >= tol)
                    {
                        signPoint = sgn(p.x - lines[i*nParam]);
                        signPoint0 = sgn(p.xprevious - lines[i*nParam]);
                        //std::cout << "signpoint1 " << signPoint << " " << signPoint0 << std::endl;
                    }
                    else
                    {
                        signPoint = sgn(p.z - p.x*lines[i*nParam+4] - lines[i*nParam+5]);
                        signPoint0 = sgn(p.zprevious - p.xprevious*lines[i*nParam+4] - lines[i*nParam+5]);
                        //std::cout << "signpoint2 " << signPoint << " " << signPoint0 << std::endl;
                    }

                    if (signPoint != signPoint0)
                    {
                        signLine1 = sgn(lines[i*nParam+1] - lines[i*nParam]*particle_slope - particle_intercept);
                        signLine2 = sgn(lines[i*nParam+3] - lines[i*nParam+2]*particle_slope - particle_intercept);
                       // std::cout << "signLines " << signLine1 << " " << signLine2 << std::endl;
                        
                        if (signLine1 != signLine2)
                        {
                            nIntersections++;
                            //std::cout << "nintersections " << nIntersections << std::endl;
                            //std::cout << fabs(p.x - p.xprevious) << tol_small << std::endl;        
                            if (fabs(p.x - p.xprevious) < tol_small)
                            {
                              //  std::cout << "vertical line" << std::cout;
                                intersectionx[nIntersections-1] = p.xprevious;
                                intersectiony[nIntersections-1] = intersectionx[nIntersections-1]*lines[i*nParam+4] +
                                                                     lines[i*nParam + 5];
                            }
                            else
                            {
                               //std::cout << "not vertical line" << std::endl;
                            //std::cout << 0.0*7.0 << " " << i << " " << nParam << " " << lines[i*nParam+4] << "  " <<tol << std::endl;
                                if (lines[i*nParam+4] >= tol)
                                {
                                    intersectionx[nIntersections-1] = lines[i*nParam];
                                }
                                else
                                {
                                    intersectionx[nIntersections-1] = (lines[i*nParam+5] - particle_intercept)/
                                                (particle_slope - lines[i*nParam+4]);
                                    //std::cout << "in this else "<< intersectionx[nIntersections -1] << std::endl;
                                }
                                intersectiony[nIntersections-1] = intersectionx[nIntersections-1]*particle_slope
                                                                               + particle_intercept;
                            }
                        }
                    }
                }

	        if (nIntersections == 0)
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
                p.x = intersectionx[0];
                //go back here and calculate p.y
                p.z = intersectiony[0];
                //std::cout << "nInt = 1 position " << intersectionx[0] << " " << intersectiony[0]  << std::endl;
            }
            else
            {
                for (int i=0; i<nIntersections; i++)
                {
                    distances[i] = (p.xprevious - intersectionx[i*nParam])*(p.xprevious - intersectionx[i*nParam]) + 
                        (p.zprevious - intersectiony[i*nParam+1])*(p.zprevious - intersectiony[i*nParam+1]);
                    if (distances[i] < minDist)
                    {
                        minDist = distances[i];
                        minDistInd = i;
                    }
                }

                p.hitWall = 1.0;
                p.x = intersectionx[minDistInd];
                p.z = intersectiony[minDistInd];
            }
    	    }
    	}
};

#endif
