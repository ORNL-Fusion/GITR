#ifndef _CFDIFFUSION_
#define _CFDIFFUSION_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particle.h"
#include <cmath>

struct crossFieldDiffusion { 

    const double dt;
	const double diffusionCoefficient;
    int nR_Bfield;
    int nZ_Bfield;
    double * BfieldGridRDevicePointer;
    double * BfieldGridZDevicePointer;
    double * BfieldRDevicePointer;
    double * BfieldZDevicePointer;
    double * BfieldTDevicePointer;

    crossFieldDiffusion(double _dt, double _diffusionCoefficient,
            int _nR_Bfield, int _nZ_Bfield,
            double * _BfieldGridRDevicePointer,double * _BfieldGridZDevicePointer,
            double * _BfieldRDevicePointer,double * _BfieldZDevicePointer,
            double * _BfieldTDevicePointer)
        : dt(_dt), diffusionCoefficient(_diffusionCoefficient),nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
       BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 

	    if(p.hitWall == 0.0)
        {
	        double perpVector[3]= {0, 0, 0};
	        double B[3] = {0.0,0.0,0.0};
            double Bmag = 0.0;
		double B_unit[3] = {0.0, 0.0, 0.0};
		double phi_random;
		double norm;
		double step;
        B[0] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
                                     BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer);
        B[2] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldZDevicePointer);
        B[1] = interp2dCombined(p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldTDevicePointer);
        Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;
#ifdef __CUDACC__
        	double r3 = curand_uniform(&p.streams[2]);
#else
        	std::uniform_real_distribution<double> dist(0.0, 1.0);
        	double r3=dist(p.streams[2]);
#endif 
		phi_random = 2*3.14159265*r3;
		perpVector[0] = cos(phi_random);
		perpVector[1] = sin(phi_random);
		perpVector[2] = (-perpVector[0]*B_unit[0] - perpVector[1]*B_unit[1])/B_unit[2];

		if (B_unit[2] == 0){
			perpVector[2] = perpVector[1];
			perpVector[1] = (-perpVector[0]*B_unit[0] - perpVector[2]*B_unit[2])/B_unit[1];
		}
		
		if ((B_unit[0] == 1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0) || (B_unit[0] == -1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0))
		{
			perpVector[2] = perpVector[0];
			perpVector[0] = 0;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==1.0 && B_unit[2] ==0.0) || (B_unit[0] == 0.0 && B_unit[1] ==-1.0 && B_unit[2] ==0.0))
		{
			perpVector[1] = 0.0;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==1.0) || (B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==-1.0))
		{
			perpVector[2] = 0;
		}
		
		norm = sqrt(perpVector[0]*perpVector[0] + perpVector[1]*perpVector[1] + perpVector[2]*perpVector[2]);
		perpVector[0] = perpVector[0]/norm;
		perpVector[1] = perpVector[1]/norm;
		perpVector[2] = perpVector[2]/norm;
		
		step = sqrt(6*diffusionCoefficient*dt);

		p.xprevious = p.xprevious + step*perpVector[0];
		p.yprevious = p.yprevious + step*perpVector[1];
		p.zprevious = p.zprevious + step*perpVector[2];
    	}
    } 
};

#endif
