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
#include "Particle.h"
#include "Boundary.h"
const double B[3] = {0.0,0.0,-2.0};

template <typename T>
CUDA_CALLABLE_MEMBER
int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

CUDA_CALLABLE_MEMBER

#ifdef __CUDACC__
void getE ( double x, double y, double z, double E[], Boundary *boundaryVector, int nLines ) {
#else
void getE ( double x, double y, double z, double E[], std::vector<Boundary> &boundaryVector, int nLines ) {
#endif
	double Emag;
	double fd = 0.8357;
	double lane = 3.64479e-04;
	double dl = 1.05058e-05;
	double pot = 60.0;
    int minIndex = 0;
    double minDistance = 1e12;
    int direction_type;
    double tol = 1e12;
    double point1_dist;
    double point2_dist;
    double perp_dist;
    double directionUnitVector[3] = {0.0,0.0,0.0};
    double vectorMagnitude;
    double max = 0.0;
    double min = 0.0;
    for (int j=0; j< nLines; j++)
    {
        if (boundaryVector[j].Z != 0.0)
        {
            point1_dist = sqrt((x - boundaryVector[j].x1)*(x - boundaryVector[j].x1) + 
                    (z - boundaryVector[j].z1)*(z - boundaryVector[j].z1));
            point2_dist = sqrt((x - boundaryVector[j].x2)*(x - boundaryVector[j].x2) + 
                                        (z - boundaryVector[j].z2)*(z - boundaryVector[j].z2));
            perp_dist = (boundaryVector[j].slope_dzdx*x - z + boundaryVector[j].intercept_z)/
                sqrt(boundaryVector[j].slope_dzdx*boundaryVector[j].slope_dzdx + 1);   

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
            if (boundaryVector[j].length*boundaryVector[j].length + perp_dist*perp_dist >=
                    max*max)
            {
                boundaryVector[j].distanceToParticle =fabs( perp_dist);
                boundaryVector[j].pointLine = 1;
            }
            else
            {
                boundaryVector[j].distanceToParticle = min;
                if (boundaryVector[j].distanceToParticle == point1_dist)
                {
                    boundaryVector[j].pointLine = 2;
                }
                else
                {
                    boundaryVector[j].pointLine = 3;
                }
            }

            if (boundaryVector[j].distanceToParticle < minDistance)
            {
                minDistance = boundaryVector[j].distanceToParticle;
                minIndex = j;
                direction_type = boundaryVector[j].pointLine;
            }
        }
        else
        {
            boundaryVector[j].distanceToParticle = tol;
        }
    if (direction_type == 1)
    {
        if (boundaryVector[minIndex].slope_dzdx == 0)
        {
            directionUnitVector[0] = 0.0;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 1.0*sgn(boundaryVector[minIndex].z1 - z);
        }
        else if (fabs(boundaryVector[minIndex].slope_dzdx)>= 0.75*tol)
        {
            
            directionUnitVector[0] = boundaryVector[minIndex].x1 - x;
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = 0.0;
        }
        else
        {
            directionUnitVector[0] = -1.0*sgn(boundaryVector[minIndex].slope_dzdx)*sgn(perp_dist);
            directionUnitVector[1] = 0.0;
            directionUnitVector[2] = -1.0*directionUnitVector[0]/(boundaryVector[minIndex].slope_dzdx);
        }
    }
    else if (direction_type == 2)
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x1 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z1 - z);
    }
    else
    {
        directionUnitVector[0] = (boundaryVector[minIndex].x2 - x);
        directionUnitVector[1] = 0.0;
        directionUnitVector[2] = (boundaryVector[minIndex].z2 - z);
    }

    vectorMagnitude = sqrt(directionUnitVector[0]*directionUnitVector[0] + directionUnitVector[1]*directionUnitVector[1]
                                + directionUnitVector[2]*directionUnitVector[2]);
    directionUnitVector[0] = directionUnitVector[0]/vectorMagnitude;
    directionUnitVector[1] = directionUnitVector[1]/vectorMagnitude;
    directionUnitVector[2] = directionUnitVector[2]/vectorMagnitude;
   
        Emag = pot*(fd/(2.0*dl)*exp(-minDistance/(2.0*dl))+ (1.0 - fd)/(lane)*exp(-minDistance/lane) );
        E[0] = Emag*directionUnitVector[0];
        E[1] = Emag*directionUnitVector[1];
        E[2] = Emag*directionUnitVector[2];

    }

}

struct move_boris { 
#ifdef __CUDACC__
    Boundary *boundaryVector;
#else
    std::vector<Boundary> &boundaryVector;
#endif
    const double span;
    const int nLines;
#ifdef __CUDACC__
    move_boris(double _span, Boundary *_boundaryVector,int _nLines) : span(_span), boundaryVector(_boundaryVector), nLines(_nLines) {}
#else
    move_boris(double _span, std::vector<Boundary> &_boundaryVector, int _nLines) : span(_span), boundaryVector(_boundaryVector), nLines(_nLines) {}
#endif    

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
	        double v_minus[3]= {0, 0, 0};
	        double v[3]= {0, 0, 0};
	        double E[3] = {0, 0, 0};
	        double B[3] = {0.0,0.0,-2.0};
	        double dt = span;
	        double Bmag = 2;
	        double q_prime = p.Z*1.60217662e-19/(p.amu*1.6737236e-27)*dt*0.5;
            double coeff = 2*q_prime/(1+(q_prime*Bmag)*(q_prime*Bmag));
            int nSteps = floor( span / dt + 0.5);
            for ( int s=0; s<nSteps; s++ ) 
            {
	          getE(p.xprevious,p.yprevious,p.zprevious,E,boundaryVector,nLines);
                   
	            v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;
                  
	            	  
	            v_minus[0] = v[0] + q_prime*E[0];
	            v_minus[1] = v[1] + q_prime*E[1];
                v_minus[2] = v[2] + q_prime*E[2];
	                                                               
                v[0] = v_minus[0] + q_prime*(v_minus[1]*B[2] - v_minus[2]*B[1]);
                v[1] = v_minus[1] + q_prime*(v_minus[2]*B[0] - v_minus[0]*B[2]);
                v[2] = v_minus[2] + q_prime*(v_minus[0]*B[1] - v_minus[1]*B[0]);
	                                    
                v[0] = v_minus[0] + coeff*(v[1]*B[2] - v[2]*B[1]);
                v[1] = v_minus[1] + coeff*(v[2]*B[0] - v[0]*B[2]);
                v[2] = v_minus[2] + coeff*(v[0]*B[1] - v[1]*B[0]);
	                                                     
	            v[0] = v[0] + q_prime*E[0];
	            v[1] = v[1] + q_prime*E[1];
	            v[2] = v[2] + q_prime*E[2];
	                                                                    
	            p.x = p.xprevious + v[0]*dt;
	            p.y = p.yprevious + v[1]*dt;
	            p.z = p.zprevious + v[2]*dt;
	                                                                    
	                             
	            						
	            p.vx = v[0];
	            p.vy = v[1];
	            p.vz = v[2];
    	    }
    	}
    } 
};

#endif
