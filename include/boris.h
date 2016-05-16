#ifndef _BORIS_
#define _BORIS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particle.h"

const double B[3] = {0.0,0.0,-2.0};

CUDA_CALLABLE_MEMBER
void getE ( double x, double y, double z, double E[] ) {

	double Emag;
	double surfaceDirection[3] = {1.7321, 0, -1.00};
	double surfaceDirection_unit[3] =  {0.866, 0, -0.50};
	double perpDistanceToSurface = ( -surfaceDirection[0]*x + z )/2.0;
	double fd = 0.8357;
	double lane = 3.64479e-04;
	double dl = 1.05058e-05;
	double pot = 60.0;
	
	Emag = pot*(fd/(2.0*dl)*exp(-perpDistanceToSurface/(2.0*dl))+ (1.0 - fd)/(lane)*exp(-perpDistanceToSurface/lane) );
	E[0] = Emag*surfaceDirection_unit[0];
	E[1] = Emag*surfaceDirection_unit[1];
	E[2] = Emag*surfaceDirection_unit[2];



}

struct move_boris { 

    const double span;

    move_boris(double _span) : span(_span) {} 

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
	        double v_minus[3]= {0, 0, 0};
	        double v[3]= {0, 0, 0};
	        double E[3] = {0, 0, 0};
	        double r[3] = {0, 0, 0};
	        double surfaceDirection[3] = {1.7321, 0, -1.00};
	        double surfaceDirection_unit[3] = {0.866, 0, -0.50};
	        double t;
	        double surface_dz_dx;
	        double * Emag;
	        double B[3] = {0.0,0.0,-2.0};
	        double perpDistanceToSurface = ( -surfaceDirection[0]*p.xprevious + p.zprevious )/2.0;	
	        double dt = 1e-9;
	        double Bmag = 2;
	        double q_prime = p.Z*1.60217662e-19/(p.amu*1.6737236e-27)*dt*0.5;
            double coeff = 2*q_prime/(1+(q_prime*Bmag)*(q_prime*Bmag));
  
            int nSteps = floor( span / dt + 0.5);

            for ( int s=0; s<nSteps; s++ ) 
            {
	            getE(p.xprevious,p.yprevious,p.zprevious,E);
	            surface_dz_dx = surfaceDirection[0];
                   
	            v[0] = p.vx;
                v[1] = p.vy;
	            v[2] = p.vz;
                  
	            r[0] = p.xprevious;
                r[1] = p.yprevious;
	            r[2] = p.zprevious;	
	            	  
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
