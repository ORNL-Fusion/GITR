#ifndef _THERMAL_
#define _THERMAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particle.h"
#include <cmath>
#include <math.h>

struct thermalForce { 

    const double dt;

    thermalForce(double _dt) : dt(_dt) {} 

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
		double MI = 1.6737236e-27;
		double alpha;
		double beta;
		double mu;
    		double amu_backgroundIons = 2.0;
		double gradTe_radial[3] = {0.0,0.0,0.0};
		double gradTe_s[3] = {0.0,0.0,0.0};
		double gradTi_radial[3] = {0.0,0.0,0.0};
		double gradTi_s[3] = {0.0,0.0,0.0};

		mu = amu_backgroundIons/(amu_backgroundIons + p.amu);
		alpha = p.Z*p.Z*0.71;
		beta =  -3*(1- mu - 5*sqrt(2.0)*pow(p.Z,2))*(1.1*pow(mu,(5/2))
			 - 0.35*pow(mu,(3/2)))/(2.6 - 2*mu+ 5.4*pow(mu,2));
		
		p.vx = p.vx + dt/(p.amu*MI)*(alpha*(gradTi_radial[0]+gradTi_s[0]) + beta*(gradTe_radial[0]+gradTe_s[0]));
		p.vy = p.vy + dt/(p.amu*MI)*(alpha*(gradTi_radial[1]+gradTi_s[1]) + beta*(gradTe_radial[1]+gradTe_s[1]));
		p.vz = p.vz + dt/(p.amu*MI)*(alpha*(gradTi_radial[2]+gradTi_s[2]) + beta*(gradTe_radial[2]+gradTe_s[2]));		
	}
    	}
     
};

#endif
