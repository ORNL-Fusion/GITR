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
    double background_amu;
    int nR_gradT;
    int nZ_gradT;
    double* gradTGridr;
    double* gradTGridz;
    double* gradTiR;
    double* gradTiZ;
    double* gradTeR;
    double* gradTeZ;

    thermalForce(double _dt, double _background_amu,int _nR_gradT, int _nZ_gradT, double* _gradTGridr, double* _gradTGridz,
            double* _gradTiR, double* _gradTiZ, double* _gradTeR, double* _gradTeZ)
        : dt(_dt), background_amu(_background_amu),nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ), gradTeR(_gradTeR), gradTeZ(_gradTeZ) {} 

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
		double MI = 1.6737236e-27;
		double alpha;
		double beta;
		double mu;
		double gradTe[3] = {0.0,0.0,0.0};
		double gradTi[3] = {0.0,0.0,0.0};
        double dv_ITG[3] = {};
        double dv_ETG[3] = {};
        gradTi[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTiR );
        gradTi[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTiZ );
        gradTi[1] = 0.0;
                    
        gradTe[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTeR );
        gradTe[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTeZ );
        gradTe[1] = 0.0;
		mu = p.amu/(background_amu + p.amu);
		alpha = p.charge*p.charge*0.71;
		beta =  -3*(1- mu - 5*sqrt(2.0)*pow(p.charge,2))*(1.1*pow(mu,(5/2))
			 - 0.35*pow(mu,(3/2)))/(2.6 - 2*mu+ 5.4*pow(mu,2));
	dv_ETG[0] = dt/(p.amu*MI)*(alpha*(gradTe[0]));
	dv_ETG[1] = dt/(p.amu*MI)*(alpha*(gradTe[1]));
	dv_ETG[2] = dt/(p.amu*MI)*(alpha*(gradTe[2]));

	dv_ITG[0] = dt/(p.amu*MI)*(beta*(gradTi[0]));
	dv_ITG[1] = dt/(p.amu*MI)*(beta*(gradTi[1]));
	dv_ITG[2] = dt/(p.amu*MI)*(beta*(gradTi[2]));
    std::cout << "mu " << mu << std::endl;
    std::cout << "alpha beta " << alpha << " " << beta << std::endl;
    std::cout << "ITG " << dv_ITG[0] << " " << dv_ITG[1] << " " << dv_ITG[2] << std::endl;
    std::cout << "ETG " << dv_ETG[0] << " " << dv_ETG[1] << " " << dv_ETG[2] << std::endl;
    std::cout << "v before thermal force " << p.vx << " " << p.vy << " " << p.vz << std::endl;
	
        p.vx = p.vx + dt/(p.amu*MI)*(alpha*(gradTe[0]) + beta*(gradTi[0]));
		p.vy = p.vy + dt/(p.amu*MI)*(alpha*(gradTe[1]) + beta*(gradTi[1]));
		p.vz = p.vz + dt/(p.amu*MI)*(alpha*(gradTe[2]) + beta*(gradTi[2]));		
	
    std::cout << "v after thermal force " << p.vx << " " << p.vy << " " << p.vz << std::endl;
        }
    	}
     
};

#endif
