#ifndef _THERMAL_
#define _THERMAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particles.h"
#include <cmath>
#include <math.h>

struct thermalForce { 
    Particles *particlesPointer;
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
            
    thermalForce(Particles *_particlesPointer, float _dt, float _background_amu,int _nR_gradT, int _nZ_gradT, float* _gradTGridr, float* _gradTGridz,
            float* _gradTiR, float* _gradTiZ, float* _gradTiT, float* _gradTeR, float* _gradTeZ,float* _gradTeT,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,
            float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,
            float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer)
        
            : particlesPointer(_particlesPointer), dt(_dt), background_amu(_background_amu),nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ),gradTiT(_gradTiT), gradTeR(_gradTeR), gradTeZ(_gradTeZ),gradTeT(_gradTeT), 
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx) const { 

	    if((particlesPointer->hitWall[indx] == 0.0) && (particlesPointer->charge[indx] > 0.0))
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
        
        interp2dVector(&gradTi[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTiR,gradTiZ, gradTiT );
                //std::cout << "Position r z" << sqrt(particlesPointer->xprevious[indx]*particlesPointer->xprevious[indx] + particlesPointer->yprevious[indx]*particlesPointer->yprevious[indx]) << " " << particlesPointer->zprevious[indx] << std::endl;            
                //std::cout << "grad Ti " << sgn(gradTi[0])*sqrt(gradTi[0]*gradTi[0] + gradTi[1]*gradTi[1]) << " " << gradTi[2] << std::endl;            
                interp2dVector(&gradTe[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_gradT,nZ_gradT,
                    gradTGridr ,gradTGridz ,gradTeR,gradTeZ, gradTeT );
		mu = particlesPointer->amu[indx]/(background_amu + particlesPointer->amu[indx]);
		alpha = particlesPointer->charge[indx]*particlesPointer->charge[indx]*0.71;
		beta =  3*(mu + 5*sqrtf(2.0)*particlesPointer->charge[indx]*particlesPointer->charge[indx]*(1.1*powf(mu,(5/2))- 0.35*powf(mu,(3/2))) - 1)/(2.6 - 2*mu+ 5.4*powf(mu,2));
	dv_ETG[0] = dt/(particlesPointer->amu[indx]*MI)*(alpha*(gradTe[0]));
	dv_ETG[1] = dt/(particlesPointer->amu[indx]*MI)*(alpha*(gradTe[1]));
	dv_ETG[2] = dt/(particlesPointer->amu[indx]*MI)*(alpha*(gradTe[2]));

	dv_ITG[0] = dt/(particlesPointer->amu[indx]*MI)*(beta*(gradTi[0]));
	dv_ITG[1] = dt/(particlesPointer->amu[indx]*MI)*(beta*(gradTi[1]));
	dv_ITG[2] = dt/(particlesPointer->amu[indx]*MI)*(beta*(gradTi[2]));
   // std::cout << "mu " << mu << std::endl;
   // std::cout << "alpha beta " << alpha << " " << beta << std::endl;
   // std::cout << "ITG " << dv_ITG[0] << " " << dv_ITG[1] << " " << dv_ITG[2] << std::endl;
   // std::cout << "ETG " << dv_ETG[0] << " " << dv_ETG[1] << " " << dv_ETG[2] << std::endl;
   // std::cout << "v before thermal force " << particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx] << std::endl;
       
   interp2dVector(&B[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);    
        Bmag = sqrtf(B[0]*B[0] + B[1]*B[1]+ B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;

        gradTiPar = B_unit[0]*gradTi[0] + B_unit[1]*gradTi[1] + B_unit[2]*gradTi[2];
        //std::cout << "gradTi Parallel " << gradTiPar << std::endl;
        //std::cout << "gradTi Parallel " << gradTi[0]<<gradTi[1]<<gradTi[2] << std::endl;
        particlesPointer->vx[indx] = particlesPointer->vx[indx] + B_unit[0]*alpha*gradTiPar;//alpha*(gradTe[0])
		particlesPointer->vy[indx] = particlesPointer->vy[indx] + B_unit[1]*alpha*gradTiPar;//alpha*(gradTe[1])
		particlesPointer->vz[indx] = particlesPointer->vz[indx] + B_unit[2]*alpha*gradTiPar;//alpha*(gradTe[2])		
        
        //particlesPointer->vx[indx] = particlesPointer->vx[indx] + (dt/(particlesPointer->amu*MI))*(  beta*(gradTi[0]));//alpha*(gradTe[0])
		//particlesPointer->vy[indx] = particlesPointer->vy[indx] + (dt/(particlesPointer->amu*MI))*(  beta*(gradTi[1]));//alpha*(gradTe[1])
		//particlesPointer->vz[indx] = particlesPointer->vz[indx] + (dt/(particlesPointer->amu*MI))*(  beta*(gradTi[2]));//alpha*(gradTe[2])		
     //   std::cout << "dv ion thermal x" << dt/(particlesPointer->amu*MI)*(  beta*(gradTi[0])) << std::endl;	
     //  std::cout << "dv ion thermal y" << dt/(particlesPointer->amu*MI)*(  beta*(gradTi[1])) << std::endl;	
     //  std::cout << "dv ion thermal z" << dt/(particlesPointer->amu*MI)*(  beta*(gradTi[2])) << std::endl;	
        //std::cout << "theta " << theta << std::endl;
       //std::cout << "v after thermal force " << particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx] << std::endl;
  
       }
    	}
     
};

#endif
