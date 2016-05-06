#ifndef _SURFACE_
#define _SURFACE_

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
double screeningLength ( double Zprojectile, double Ztarget ) {
	double bohrRadius = 5.29177e-11;
	double screenLength;

	screenLength = 0.885341*bohrRadius*pow(pow(Zprojectile,(2/3)) + pow(Ztarget,(2/3)),(-1/2));

	return screenLength;
}

CUDA_CALLABLE_MEMBER
double stoppingPower (Particle p, double Mtarget, double Ztarget, double screenLength) {
	        double E0;
                double Q = 1.60217662e-19;
		double ke2 = 14.4e-10;
		double reducedEnergy;
	double stoppingPower;

	E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
	reducedEnergy = E0*(Mtarget/(p.amu+Mtarget))*(screenLength/(p.Z*Ztarget*ke2));
	stoppingPower = 0.5*log(1.0 + 1.2288*reducedEnergy)/(reducedEnergy + 0.1728*sqrt(reducedEnergy) + 0.008*pow(reducedEnergy, 0.1504));

	return stoppingPower;	
}

struct erosion { 

    const double dt;

    erosion(double _dt) : dt(_dt) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 
	double screenLength;
	double stopPower;
	double q = 18.6006;
	double lambda = 2.2697;
	double mu = 3.1273;
	double Eth = 24.9885;
	double Y0;
	double Ztarget = 74.0;
	double Mtarget = 183.84;
	double term;
	double E0;

	screenLength = screeningLength(p.Z, Ztarget);
	stopPower = stoppingPower(p, Mtarget, Ztarget, screenLength); 
	E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
	term = pow((E0/Eth - 1),mu);
	Y0 = q*stopPower*term/(lambda + term);
    	}
     
};

struct reflection {

    const double dt;

    reflection(double _dt) : dt(_dt) {}

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(Particle &p) const {
            if(p.hitWall == 1.0)
        	{
			double reducedEnergyMultiplier = 5e-7;//for W on W
        		double E0;
			double reducedEnergy;
			double a1 = -3.685;
			double a2 = 0.0278;
			double a3 = 7.825e-5;
			double a4 = -1.109;
			double Rn;
		
		        E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
			reducedEnergy = E0*reducedEnergyMultiplier;

			Rn = exp(a1*pow(reducedEnergy,a2))/(1 + exp(a3*pow(reducedEnergy,a4)));

#ifdef __CUDACC__
                curandState tmpState;
                double r7 = curand_uniform(&p.streams[6]);
#else
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                double r7 = dist(p.streams[6]);
#endif

#endif
		        if(r7 <= Rn)
        		{
				double a1_e = -7.168;
                        	double a2_e = 0.01685;
                        	double a3_e = 7.005e-5;
                        	double a4_e = -1.343;
				double Re = exp(a1_e*pow(reducedEnergy,a2_e))/(1 + exp(a3_e*pow(reducedEnergy,a4_e))); 
				double launch_unitVector[3] = {-0.866, 0, 0.50};
				double launchEnergy = E0*Re/Rn;
                		p.Z = 0.0;
				p.vx = launchEnergy*launch_unitVector[0];
				p.vy = launchEnergy*launch_unitVector[1];
				p.vz = launchEnergy*launch_unitVector[2];
				p.hitWall = 0.0;
        		}
		}
	}
};
