#ifndef _THERMAL_
#define _THERMAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particles.h"
#include <cmath>

struct thermalForce { 

    Particles *p;
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
	    float dv_ITGx=0.0;
	    float dv_ITGy=0.0;
	    float dv_ITGz=0.0;
	    float dv_ETGx=0.0;
	    float dv_ETGy=0.0;
	    float dv_ETGz=0.0;
            
    thermalForce(Particles *_p,float _dt, float _background_amu,int _nR_gradT, int _nZ_gradT, float* _gradTGridr, float* _gradTGridz,
            float* _gradTiR, float* _gradTiZ, float* _gradTiT, float* _gradTeR, float* _gradTeZ,float* _gradTeT,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,
            float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,
            float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer)
        
            : p(_p), dt(_dt), background_amu(_background_amu),nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ),gradTiT(_gradTiT), gradTeR(_gradTeR), gradTeZ(_gradTeZ),gradTeT(_gradTeT), 
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx)  { 
    if ((p->hitWall[indx] == 0.0) && (p->charge[indx] > 0.0)) {
      float MI = 1.6737236e-27;
      float alpha;
      float beta;
      float mu;
      float gradTe[3] = {0.0, 0.0, 0.0};
      float gradTi[3] = {0.0, 0.0, 0.0};
      float B[3] = {0.0, 0.0, 0.0};
      float B_unit[3] = {0.0, 0.0, 0.0};
      float Bmag = 0.0;
      float gradTiPar = 0.0;
      float dv_ITG[3] = {};
      float dv_ETG[3] = {};
      float vNorm = 0.0;
      float vNorm2 = 0.0;
      // std:cout << " grad Ti interp " << std::endl;
      interp2dVector(&gradTi[0], p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT,
                     gradTGridr, gradTGridz, gradTiR, gradTiZ, gradTiT);
      //std::cout << "Position r z" << sqrt(p->xprevious*p->xprevious + p->yprevious*p->yprevious) << " " << p->zprevious << std::endl;
      //std::cout << "grad Ti " << std::copysign(1.0,gradTi[0])*sqrt(gradTi[0]*gradTi[0] + gradTi[1]*gradTi[1]) << " " << gradTi[2] << std::endl;
      interp2dVector(&gradTe[0], p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT,
                     gradTGridr, gradTGridz, gradTeR, gradTeZ, gradTeT);
      mu = p->amu[indx] / (background_amu + p->amu[indx]);
      alpha = p->charge[indx] * p->charge[indx] * 0.71;
      beta = 3 * (mu + 5 * std::sqrt(2.0) * p->charge[indx] * p->charge[indx] * (1.1 * std::pow(mu, (5 / 2)) - 0.35 * std::pow(mu, (3 / 2))) - 1) / (2.6 - 2 * mu + 5.4 * std::pow(mu, 2));
       
       interp2dVector(&B[0],p->xprevious[indx],p->yprevious[indx],p->zprevious[indx],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);    
        Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1]+ B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;

	dv_ETG[0] = 1.602e-19*dt/(p->amu[indx]*MI)*(alpha*(gradTe[0]));
	dv_ETG[1] = 1.602e-19*dt/(p->amu[indx]*MI)*(alpha*(gradTe[1]));
	dv_ETG[2] = 1.602e-19*dt/(p->amu[indx]*MI)*(alpha*(gradTe[2]));

	dv_ITG[0] = 1.602e-19*dt/(p->amu[indx]*MI)*(beta*(gradTi[0]))*B_unit[0];
	dv_ITG[1] = 1.602e-19*dt/(p->amu[indx]*MI)*(beta*(gradTi[1]))*B_unit[1];
	dv_ITG[2] = 1.602e-19*dt/(p->amu[indx]*MI)*(beta*(gradTi[2]))*B_unit[2];
	dv_ITGx = dv_ITG[0];
	dv_ITGy = dv_ITG[1];
	dv_ITGz = dv_ITG[2];

    //std::cout << "mu " << mu << std::endl;
    //std::cout << "alpha beta " << alpha << " " << beta << std::endl;
    //std::cout << "ITG " << dv_ITG[0] << " " << dv_ITG[1] << " " << dv_ITG[2] << std::endl;
    //std::cout << "gradTi " << gradTi[0] << " " << gradTi[1] << " " << gradTi[2] << std::endl;
    //std::cout << "ETG " << dv_ETG[0] << " " << dv_ETG[1] << " " << dv_ETG[2] << std::endl;
    //std::cout << "v before thermal force " << p->vx[indx] << " " << p->vy[indx] << " " << p->vz[indx] << std::endl;
    /*
    float theta = atan2(p->yprevious,p->xprevious);
    float Ar = -1;
    float At = 0.0;
    float Az = 1;
    gradTi[0] = cos(theta)*Ar - sin(theta)*At;
    gradTi[1] = sin(theta)*Ar + cos(theta)*At;
    gradTi[2] = Az;
    */
    float vx = p->vx[indx];
    float vy = p->vy[indx];
    float vz = p->vz[indx];
        vNorm = std::sqrt(vx*vx + vy*vy + vz*vz);
    p->vD[indx] = dv_ITG[2];    
	//std::cout << "gradTi Parallel " << gradTiPar << std::endl;
        //std::cout << "gradTi Parallel " << gradTi[0]<<gradTi[1]<<gradTi[2] << std::endl;
        //p->vx[indx] = p->vx[indx] +dv_ITG[0];//alpha*(gradTe[0])   
	//p->vy[indx] = p->vy[indx] +dv_ITG[1];//alpha*(gradTe[1])
	//p->vz[indx] = p->vz[indx] +dv_ITG[2];//alpha*(gradTe[2])		
        //vNorm2 = sqrt(p->vx[indx]*p->vx[indx] + p->vy[indx]*p->vy[indx] + p->vz[indx]*p->vz[indx]);
		//SFT
        float k1 = dv_ITG[2] - dt*p->nu_s[indx]
                    *(dv_ITG[2]);
        p->vx[indx] = vx + dv_ITG[0];///velocityCollisionsNorm;   	
		p->vy[indx] = vy + dv_ITG[1];///velocityCollisionsNorm;   	
		p->vz[indx] = vz + dv_ITG[2];// - dt*p->nu_s[indx]
                         //*(0.5*k1);///velocityCollisionsNorm;   	
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
