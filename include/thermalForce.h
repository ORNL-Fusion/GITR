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

    const float dt;
    float background_amu;
    int nR_gradT;
    int nZ_gradT;
    float* gradTGridr;
    float* gradTGridz;
    float* gradTiR;
    float* gradTiZ;
    float* gradTiT;
    float* gradTeR;
    float* gradTeZ;
    float* gradTeT;
            int nR_Bfield;
            int nZ_Bfield;
            float * BfieldGridRDevicePointer;
            float * BfieldGridZDevicePointer;
            float * BfieldRDevicePointer;
            float * BfieldZDevicePointer;
            float * BfieldTDevicePointer;
            
    thermalForce(float _dt, float _background_amu,int _nR_gradT, int _nZ_gradT, float* _gradTGridr, float* _gradTGridz,
            float* _gradTiR, float* _gradTiZ, float* _gradTiT, float* _gradTeR, float* _gradTeZ,float* _gradTeT,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,
            float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,
            float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer)
        
            : dt(_dt), background_amu(_background_amu),nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ),gradTiT(_gradTiT), gradTeR(_gradTeR), gradTeZ(_gradTeZ),gradTeT(_gradTeT), 
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {}

CUDA_CALLABLE_MEMBER    
void operator()(Particle &p) const { 

	    if((p.hitWall == 0.0) && (p.charge > 0.0))
        {
		float MI = 1.6737236e-27;
		float alpha;
		float beta;
		float mu;
		float gradTe[3] = {0.0,0.0,0.0};
		float gradTi[3] = {0.0,0.0,0.0};
		float B[3] = {0.0,0.0,0.0};
		float B_unit[3] = {0.0,0.0,0.0};
        float Bmag = 0.0;
        float gradTiPar = 0.0;
        float dv_ITG[3] = {};
        float dv_ETG[3] = {};
               // std:cout << " grad Ti interp " << std::endl;
                interp2dVector(&gradTi[0],p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTiR,gradTiZ, gradTiT );
                //std::cout << "Position r z" << sqrt(p.xprevious*p.xprevious + p.yprevious*p.yprevious) << " " << p.zprevious << std::endl;            
                //std::cout << "grad Ti " << sgn(gradTi[0])*sqrt(gradTi[0]*gradTi[0] + gradTi[1]*gradTi[1]) << " " << gradTi[2] << std::endl;            
                interp2dVector(&gradTe[0],p.xprevious,p.yprevious,p.zprevious,nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTeR,gradTeZ, gradTeT );
		mu = p.amu/(background_amu + p.amu);
		alpha = p.charge*p.charge*0.71;
		beta =  3*(mu + 5*sqrt(2.0)*p.charge*p.charge*(1.1*pow(mu,(5/2))- 0.35*pow(mu,(3/2))) - 1)/(2.6 - 2*mu+ 5.4*pow(mu,2));
	dv_ETG[0] = dt/(p.amu*MI)*(alpha*(gradTe[0]));
	dv_ETG[1] = dt/(p.amu*MI)*(alpha*(gradTe[1]));
	dv_ETG[2] = dt/(p.amu*MI)*(alpha*(gradTe[2]));

	dv_ITG[0] = dt/(p.amu*MI)*(beta*(gradTi[0]));
	dv_ITG[1] = dt/(p.amu*MI)*(beta*(gradTi[1]));
	dv_ITG[2] = dt/(p.amu*MI)*(beta*(gradTi[2]));
   // std::cout << "mu " << mu << std::endl;
   // std::cout << "alpha beta " << alpha << " " << beta << std::endl;
   // std::cout << "ITG " << dv_ITG[0] << " " << dv_ITG[1] << " " << dv_ITG[2] << std::endl;
   // std::cout << "ETG " << dv_ETG[0] << " " << dv_ETG[1] << " " << dv_ETG[2] << std::endl;
   // std::cout << "v before thermal force " << p.vx << " " << p.vy << " " << p.vz << std::endl;
    /*
    float theta = atan2(p.yprevious,p.xprevious);
    float Ar = -1;
    float At = 0.0;
    float Az = 1;
    gradTi[0] = cos(theta)*Ar - sin(theta)*At;
    gradTi[1] = sin(theta)*Ar + cos(theta)*At;
    gradTi[2] = Az;
    */
       interp2dVector(&B[0],p.xprevious,p.yprevious,p.zprevious,nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);    
        Bmag = sqrt(B[0]*B[0] + B[1]*B[1]+ B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;

        gradTiPar = B_unit[0]*gradTi[0] + B_unit[1]*gradTi[1] + B_unit[2]*gradTi[2];
        //std::cout << "gradTi Parallel " << gradTiPar << std::endl;
        //std::cout << "gradTi Parallel " << gradTi[0]<<gradTi[1]<<gradTi[2] << std::endl;
        p.vx = p.vx + B_unit[0]*alpha*gradTiPar;//alpha*(gradTe[0])
		p.vy = p.vy + B_unit[1]*alpha*gradTiPar;//alpha*(gradTe[1])
		p.vz = p.vz + B_unit[2]*alpha*gradTiPar;//alpha*(gradTe[2])		
        //p.vx = p.vx + (dt/(p.amu*MI))*(  beta*(gradTi[0]));//alpha*(gradTe[0])
		//p.vy = p.vy + (dt/(p.amu*MI))*(  beta*(gradTi[1]));//alpha*(gradTe[1])
		//p.vz = p.vz + (dt/(p.amu*MI))*(  beta*(gradTi[2]));//alpha*(gradTe[2])		
     //   std::cout << "dv ion thermal x" << dt/(p.amu*MI)*(  beta*(gradTi[0])) << std::endl;	
     //  std::cout << "dv ion thermal y" << dt/(p.amu*MI)*(  beta*(gradTi[1])) << std::endl;	
     //  std::cout << "dv ion thermal z" << dt/(p.amu*MI)*(  beta*(gradTi[2])) << std::endl;	
        //std::cout << "theta " << theta << std::endl;
       //std::cout << "v after thermal force " << p.vx << " " << p.vy << " " << p.vz << std::endl;
        }
    	}
     
};

#endif
