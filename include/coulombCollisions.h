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
			 	double& nu_energy, Particle& p, 
                
    int nR_flowV,
    int nZ_flowV,
    double* flowVGridr,
    double* flowVGridz,
    double* flowVr,
    double* flowVz,
    double* flowVt,
    int nR_Dens,
    int nZ_Dens,
    double* DensGridr,
    double* DensGridz,
    double* ni,
    int nR_Temp,
    int nZ_Temp,
    double* TempGridr,
    double* TempGridz,
    double* ti, double background_Z, double background_amu
                ) {
        double Q = 1.60217662e-19;
        double EPS0 = 8.854187e-12;
	double pi = 3.14159265;
	
        double Temp_eV = interp2dCombined(p.x,p.y,p.z,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti);
            double density = interp2dCombined(p.x,p.y,p.z,nR_Dens,nZ_Dens,DensGridr,DensGridz,ni);
	double flowVelocity[3]= {0, 0, 0};
	double relativeVelocity[3] = {0.0, 0.0, 0.0};
	double velocityNorm;
	double lam_d;
	double lam;
	double gam;
	double a;
	double x;
	double psi_prime;
	double psi_psiprime;
	double psi;
	double nu_0;
    flowVelocity[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVr);
    flowVelocity[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVz);
    flowVelocity[1] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVt);
	relativeVelocity[0] = p.vx - flowVelocity[0];
	relativeVelocity[1] = p.vy - flowVelocity[1];
	relativeVelocity[2] = p.vz - flowVelocity[2];
	velocityNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);                
		//for(int i=1; i < nSpecies; i++)
		//{
			lam_d = sqrt(EPS0*Temp_eV/(density*pow(background_Z,2)*Q));//only one q in order to convert to J
                	lam = 4.0*pi*density*pow(lam_d,3);
                	gam = pow(Q,4)*pow(p.Z,2)*pow(background_Z,2)*log(lam)/(p.amu*p.amu*4*pi*EPS0*EPS0);
                
                	a = background_amu/(2*Temp_eV*Q);// %q is just to convert units - no z needed
                
                	x = pow(velocityNorm,2)*a;
                	psi_prime = 2*sqrt(x/pi)*exp(-x);
                	psi_psiprime = erf(sqrt(x));
                	psi = psi_psiprime - psi_prime;
                	nu_0 = gam*density/pow(velocityNorm,3);
                	nu_friction = nu_friction -(1+p.amu/background_amu)*psi*nu_0;
                	nu_deflection = nu_deflection + 2*(psi_psiprime - psi/(2*x))*nu_0;
                	nu_parallel = nu_parallel + psi/x*nu_0;
                	nu_energy = nu_energy+2*(p.amu/background_amu*psi - psi_prime)*nu_0;
	//	}
}

CUDA_CALLABLE_MEMBER
void getSlowDownDirections (double parallel_direction[], double perp_direction1[], double perp_direction2[], Particle& p,
    int nR_flowV,
    int nZ_flowV,
    double* flowVGridr,
    double* flowVGridz,
    double* flowVr,
    double* flowVz,
    double* flowVt,
    
                        int nR_Bfield, int nZ_Bfield,
                        double* BfieldGridR ,double* BfieldGridZ ,
                        double* BfieldR ,double* BfieldZ ,
                 double* BfieldT 
    ) {
	        double flowVelocity[3]= {0, 0, 0.0};
                double relativeVelocity[3] = {0.0, 0.0, 0.0};
                double B[3] = {};
                double Bmag = 0.0;
		double B_unit[3] = {0.0, 0.0, -1.0};
		double velocityRelativeNorm;
		double s1;
		double s2;
        B[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                                                     BfieldGridR ,BfieldGridZ ,BfieldR );
        B[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                               BfieldGridR ,BfieldGridZ ,BfieldZ );
        B[1] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                               BfieldGridR ,BfieldGridZ ,BfieldT );
        Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;
        flowVelocity[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVr);
    flowVelocity[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVz);
    flowVelocity[1] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
        flowVGridr,flowVGridz,flowVt);

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
    int nR_flowV;
    int nZ_flowV;
    double* flowVGridr;
    double* flowVGridz;
    double* flowVr;
    double* flowVz;
    double* flowVt;
    int nR_Dens;
    int nZ_Dens;
    double* DensGridr;
    double* DensGridz;
    double* ni;
    int nR_Temp;
    int nZ_Temp;
    double* TempGridr;
    double* TempGridz;
    double* ti;
    double background_Z;
    double background_amu;
    int nR_Bfield;
    int nZ_Bfield;
    double * BfieldGridR;
    double * BfieldGridZ;
    double * BfieldR;
    double * BfieldZ;
    double * BfieldT;
    coulombCollisions(double _dt, int _nR_flowV, int _nZ_flowV,    double* _flowVGridr,
                double* _flowVGridz,double* _flowVr,
                        double* _flowVz,double* _flowVt,
                        int _nR_Dens,int _nZ_Dens,double* _DensGridr,
                            double* _DensGridz,double* _ni,int _nR_Temp, int _nZ_Temp,
                        double* _TempGridr, double* _TempGridz,double* _ti,
                        double _background_Z, double _background_amu,
                        int _nR_Bfield, int _nZ_Bfield,
                        double * _BfieldGridR ,double * _BfieldGridZ ,
                        double * _BfieldR ,double * _BfieldZ ,
                 double * _BfieldT )
        : dt(_dt), nR_flowV(_nR_flowV), nZ_flowV(_nZ_flowV), flowVGridr(_flowVGridr),
   flowVGridz(_flowVGridz), flowVr(_flowVr),flowVz(_flowVz), flowVt(_flowVt),
   nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),ni(_ni),
           nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz),
           ti(_ti),background_Z(_background_Z), background_amu(_background_amu),
   nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridR(_BfieldGridR), 
    BfieldGridZ(_BfieldGridZ),BfieldR(_BfieldR), BfieldZ(_BfieldZ), BfieldT(_BfieldT) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
		double nu_friction = 0.0;
		double nu_deflection = 0.0;
		double nu_parallel = 0.0;
		double nu_energy = 0.0;
		double flowVelocity[3]= {0, 0, 0.0};
		double relativeVelocity[3] = {0.0, 0.0, 0.0};
		double velocityCollisions[3];	
		double velocityRelativeNorm;	
		double parallel_direction[3];
		double perp_direction1[3];
		double perp_direction2[3];
		double parallel_contribution;
		double dv_perp1[3];
		double dv_perp2[3];
	
        flowVelocity[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
                     flowVGridr,flowVGridz,flowVr);
        flowVelocity[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
                     flowVGridr,flowVGridz,flowVz);
        flowVelocity[1] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_flowV,nZ_flowV,
                     flowVGridr,flowVGridz,flowVt);
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

		getSlowDownFrequencies (nu_friction,nu_deflection,nu_parallel, nu_energy,p ,
                  nR_flowV,  nZ_flowV, flowVGridr,
        flowVGridz,flowVr,
        flowVz,flowVt,
         nR_Dens, nZ_Dens,DensGridr,
        DensGridz,ni, nR_Temp,  nZ_Temp,
        TempGridr, TempGridz,ti, background_Z, background_amu);
		getSlowDownDirections(parallel_direction, perp_direction1,perp_direction2,p,
                  nR_flowV,  nZ_flowV, flowVGridr,
        flowVGridz,flowVr,
        flowVz,flowVt,
                
    nR_Bfield,
    nZ_Bfield,
    BfieldGridR,
    BfieldGridZ,
    BfieldR,
    BfieldZ,
    BfieldT
                );
		
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
