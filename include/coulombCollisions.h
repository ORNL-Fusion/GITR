#ifndef _COULOMB_
#define _COULOMB_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particle.h"
#include <cmath>
#include <math.h>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#include <stdlib.h>
#endif

CUDA_CALLABLE_MEMBER
void getSlowDownFrequencies ( double& nu_friction, double& nu_deflection, double& nu_parallel,
			 	double& nu_energy, Particle& p ) {
        double Q = 1.60217662e-19;
        double EPS0 = 8.854187e-12;
	double pi = 3.14159265;
	int nSpecies = 2;
	double density = 1e19;
        double Temp_eV = 20;
	double flowVelocity[3]= {0, 0, -1000.0};
	double relativeVelocity[3] = {0.0, 0.0, 0.0};
	double velocityNorm;
	double backgroundZ[2] = {-1.0, 1.0};
	double backgroundAMU[2] = {0.000548597,2.0};
	double lam_d;
	double lam;
	double gam;
	double a;
	double x;
	double psi_prime;
	double psi_psiprime;
	double psi;
	double nu_0;

	relativeVelocity[0] = p.vx - flowVelocity[0];
	relativeVelocity[1] = p.vy - flowVelocity[1];
	relativeVelocity[2] = p.vz - flowVelocity[2];
	velocityNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);                
		for(int i=1; i < nSpecies; i++)
		{
			lam_d = sqrt(EPS0*Temp_eV/(density*pow(backgroundZ[i],2)*Q));//only one q in order to convert to J
                	lam = 4.0*pi*density*pow(lam_d,3);
                	gam = pow(Q,4)*pow(p.Z,2)*pow(backgroundZ[i],2)*log(lam)/(p.amu*p.amu*4*pi*EPS0*EPS0);
                
                	a = backgroundAMU[i]/(2*Temp_eV*Q);// %q is just to convert units - no z needed
                
                	x = pow(velocityNorm,2)*a;
                	psi_prime = 2*sqrt(x/pi)*exp(-x);
                	psi_psiprime = erf(sqrt(x));
                	psi = psi_psiprime - psi_prime;
                	nu_0 = gam*density/pow(velocityNorm,3);
                	nu_friction = nu_friction -(1+p.amu/backgroundAMU[i])*psi*nu_0;
                	nu_deflection = nu_deflection + 2*(psi_psiprime - psi/(2*x))*nu_0;
                	nu_parallel = nu_parallel + psi/x*nu_0;
                	nu_energy = nu_energy+2*(p.amu/backgroundAMU[i]*psi - psi_prime)*nu_0;
		}
}

CUDA_CALLABLE_MEMBER
void getSlowDownDirections (double parallel_direction[], double perp_direction1[], double perp_direction2[], Particle& p) {
	        double flowVelocity[3]= {0, 0, -1000.0};
                double relativeVelocity[3] = {0.0, 0.0, 0.0};
		double B_unit[3] = {0.0, 0.0, -1.0};
		double velocityRelativeNorm;
		double s1;
		double s2;

		relativeVelocity[0] = p.vx - flowVelocity[0];
                relativeVelocity[1] = p.vy - flowVelocity[1];
                relativeVelocity[2] = p.vz - flowVelocity[2];
                velocityRelativeNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);

		parallel_direction[0] = relativeVelocity[0]/velocityRelativeNorm;
		parallel_direction[1] = relativeVelocity[1]/velocityRelativeNorm;
		parallel_direction[2] = relativeVelocity[2]/velocityRelativeNorm;

		s1 = parallel_direction[0]*B_unit[0]+parallel_direction[1]*B_unit[1]+parallel_direction[2]*B_unit[2];
            	s2 = sqrt(1-s1*s1);
            
            	perp_direction1[0] = 1/s2*(s1*parallel_direction[0] - B_unit[0]);
		perp_direction1[1] = 1/s2*(s1*parallel_direction[1] - B_unit[1]);
		perp_direction1[2] = 1/s2*(s1*parallel_direction[2] - B_unit[2]);
           
                perp_direction2[0] = 1/s2*(parallel_direction[1]*B_unit[2] - parallel_direction[2]*B_unit[1]);
                perp_direction2[1] = 1/s2*(parallel_direction[2]*B_unit[0] - parallel_direction[0]*B_unit[2]);
                perp_direction2[2] = 1/s2*(parallel_direction[0]*B_unit[1] - parallel_direction[1]*B_unit[0]); 
            if (s2 == 0)
		{
                perp_direction2[0] = parallel_direction[2];
		perp_direction2[1] = parallel_direction[0];
		perp_direction2[2] =  parallel_direction[1];

                s1 = parallel_direction[0]*perp_direction2[0]+parallel_direction[1]*perp_direction2[1]+parallel_direction[2]*perp_direction2[2];
                s2 = sqrt(1-s1*s1);
		perp_direction1[0] = -1/s2*(parallel_direction[1]*perp_direction2[2] - parallel_direction[2]*perp_direction2[1]);
                perp_direction1[1] = -1/s2*(parallel_direction[2]*perp_direction2[0] - parallel_direction[0]*perp_direction2[2]);
                perp_direction1[2] = -1/s2*(parallel_direction[0]*perp_direction2[1] - parallel_direction[1]*perp_direction2[0]);
            }
	
}

struct coulombCollisions { 

    const double dt;

    coulombCollisions(double _dt) : dt(_dt) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
		double nu_friction = 0.0;
		double nu_deflection = 0.0;
		double nu_parallel = 0.0;
		double nu_energy = 0.0;
		double flowVelocity[3]= {0, 0, -1000.0};
		double relativeVelocity[3] = {0.0, 0.0, 0.0};
		double velocityCollisions[3];	
		double velocityRelativeNorm;	
		double parallel_direction[3];
		double perp_direction1[3];
		double perp_direction2[3];
		double parallel_contribution;
		double dv_perp1[3];
		double dv_perp2[3];
	
	        relativeVelocity[0] = p.vx - flowVelocity[0];
        	relativeVelocity[1] = p.vy - flowVelocity[1];
        	relativeVelocity[2] = p.vz - flowVelocity[2];
        	velocityRelativeNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);
#ifdef __CUDACC__
        	int plus_minus1 = floor(curand_uniform(&p.streams[3]) + 0.5)*2 -1;
		int plus_minus2 = floor(curand_uniform(&p.streams[4]) + 0.5)*2 -1;
		int plus_minus3 = floor(curand_uniform(&p.streams[5]) + 0.5)*2 -1;
#else
	        std::uniform_real_distribution<double> dist(0.0, 1.0);
        	int plus_minus1 = floor(dist(p.streams[3]) + 0.5)*2 - 1;
		int plus_minus2 = floor(dist(p.streams[4]) + 0.5)*2 - 1;
		int plus_minus3 = floor(dist(p.streams[5]) + 0.5)*2 - 1;
#endif

		getSlowDownFrequencies (nu_friction,nu_deflection,nu_parallel, nu_energy,p );
		getSlowDownDirections(parallel_direction, perp_direction1,perp_direction2,p);
		
		parallel_contribution = (1.0+nu_friction*dt + plus_minus1*sqrt(nu_parallel*dt));
		dv_perp1[0] = perp_direction1[0]*plus_minus2*sqrt(nu_deflection/2*dt);
		dv_perp1[1] = perp_direction1[1]*plus_minus2*sqrt(nu_deflection/2*dt);
		dv_perp1[2] = perp_direction1[2]*plus_minus2*sqrt(nu_deflection/2*dt);
                dv_perp2[0] = perp_direction2[0]*plus_minus3*sqrt(nu_deflection/2*dt);
                dv_perp2[1] = perp_direction2[1]*plus_minus3*sqrt(nu_deflection/2*dt);
                dv_perp2[2] = perp_direction2[2]*plus_minus3*sqrt(nu_deflection/2*dt);

		velocityCollisions[0] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[0]*parallel_contribution
 					+ dv_perp1[0] + dv_perp2[0]);
		velocityCollisions[1] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[1]*parallel_contribution                                         + dv_perp1[1] + dv_perp2[1]);
		velocityCollisions[2] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[2]*parallel_contribution                                         + dv_perp1[2] + dv_perp2[2]);

		p.vx = velocityCollisions[0] + flowVelocity[0]; 
		p.vy = velocityCollisions[1] + flowVelocity[1];
		p.vz = velocityCollisions[2] + flowVelocity[2];   	
    
	}
    	}
     
};

#endif
