#ifndef _THERMAL_
#define _THERMAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particles.h"
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct thermalForce { 
    Flags* flags;
    Particles *p;
    const gitr_precision dt;
    gitr_precision background_amu;
    int nR_gradT;
    int nZ_gradT;
    gitr_precision* gradTGridr;
    gitr_precision* gradTGridz;
    gitr_precision* gradTiR;
    gitr_precision* gradTiZ;
    gitr_precision* gradTiT;
    gitr_precision* gradTeR;
    gitr_precision* gradTeZ;
    gitr_precision* gradTeT;
            int nR_Bfield;
            int nZ_Bfield;
            gitr_precision * BfieldGridRDevicePointer;
            gitr_precision * BfieldGridZDevicePointer;
            gitr_precision * BfieldRDevicePointer;
            gitr_precision * BfieldZDevicePointer;
            gitr_precision * BfieldTDevicePointer;
	    gitr_precision dv_ITGx=0.0;
	    gitr_precision dv_ITGy=0.0;
	    gitr_precision dv_ITGz=0.0;
	    gitr_precision dv_ETGx=0.0;
	    gitr_precision dv_ETGy=0.0;
	    gitr_precision dv_ETGz=0.0;
      int cylsymm;
            
    thermalForce(Flags* _flags,Particles *_p,gitr_precision _dt, gitr_precision _background_amu,int _nR_gradT, int _nZ_gradT, gitr_precision* _gradTGridr, gitr_precision* _gradTGridz,
            gitr_precision* _gradTiR, gitr_precision* _gradTiZ, gitr_precision* _gradTiT, gitr_precision* _gradTeR, gitr_precision* _gradTeZ,gitr_precision* _gradTeT,
            int _nR_Bfield, int _nZ_Bfield,
            gitr_precision * _BfieldGridRDevicePointer,
            gitr_precision * _BfieldGridZDevicePointer,
            gitr_precision * _BfieldRDevicePointer,
            gitr_precision * _BfieldZDevicePointer,
            gitr_precision * _BfieldTDevicePointer,
            int cylsymm_ )
        
            : flags(_flags),p(_p), dt(_dt), background_amu(_background_amu),nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ),gradTiT(_gradTiT), gradTeR(_gradTeR), gradTeZ(_gradTeZ),gradTeT(_gradTeT), 
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
    BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer), cylsymm( cylsymm_ ) {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx)  { 
    if ((p->hitWall[indx] == 0.0) && (p->charge[indx] > 0.0)) {
      gitr_precision MI = 1.6737236e-27;
      gitr_precision alpha;
      gitr_precision beta;
      gitr_precision mu;
      gitr_precision gradTe[3] = {0.0, 0.0, 0.0};
      gitr_precision gradTi[3] = {0.0, 0.0, 0.0};
      gitr_precision B[3] = {0.0, 0.0, 0.0};
      gitr_precision B_unit[3] = {0.0, 0.0, 0.0};
      gitr_precision Bmag = 0.0;
      gitr_precision gradTiPar = 0.0;
      gitr_precision dv_ITG[3] = {};
      gitr_precision dv_ETG[3] = {};
      gitr_precision vNorm = 0.0;
      gitr_precision vNorm2 = 0.0;
      gitr_precision dt_step = dt;
                if (flags->USE_ADAPTIVE_DT) {
	          if(p->advance[indx])
		  {
	            dt_step = p->dt[indx];
		  }
		  else
		  {
	            dt_step = 0.0;
		  }
                }
      // std:cout << " grad Ti interp " << std::endl;
      interp2dVector(&gradTi[0], p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT,
                     gradTGridr, gradTGridz, gradTiR, gradTiZ, gradTiT, cylsymm );
      //std::cout << "Position r z" << sqrt(p->xprevious*p->xprevious + p->yprevious*p->yprevious) << " " << p->zprevious << std::endl;
      //std::cout << "grad Ti " << std::copysign(1.0,gradTi[0])*sqrt(gradTi[0]*gradTi[0] + gradTi[1]*gradTi[1]) << " " << gradTi[2] << std::endl;
      interp2dVector(&gradTe[0], p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT,
                     gradTGridr, gradTGridz, gradTeR, gradTeZ, gradTeT, cylsymm );
      mu = p->amu[indx] / (background_amu + p->amu[indx]);
      alpha = p->charge[indx] * p->charge[indx] * 0.71;
      beta = 3 * (mu + 5 * std::sqrt(2.0) * p->charge[indx] * p->charge[indx] * (1.1 * std::pow(mu, (5 / 2)) - 0.35 * std::pow(mu, (3 / 2))) - 1) / (2.6 - 2 * mu + 5.4 * std::pow(mu, 2));
       
       interp2dVector(&B[0],p->xprevious[indx],p->yprevious[indx],p->zprevious[indx],nR_Bfield,nZ_Bfield,
             BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer, cylsymm );    
        Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1]+ B[2]*B[2]);
        /* Captain! These are not checked for zero! Check them for zero! */
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;

	dv_ETG[0] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(alpha*(gradTe[0]))*B_unit[0];
	dv_ETG[1] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(alpha*(gradTe[1]))*B_unit[1];
	dv_ETG[2] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(alpha*(gradTe[2]))*B_unit[2];

	dv_ETGx = dv_ETG[0];
	dv_ETGy = dv_ETG[1];
	dv_ETGz = dv_ETG[2];

	dv_ITG[0] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(beta*(gradTi[0]))*B_unit[0];
	dv_ITG[1] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(beta*(gradTi[1]))*B_unit[1];
	dv_ITG[2] = 1.602e-19*dt_step/(p->amu[indx]*MI)*(beta*(gradTi[2]))*B_unit[2];

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
    gitr_precision theta = atan2(p->yprevious,p->xprevious);
    gitr_precision Ar = -1;
    gitr_precision At = 0.0;
    gitr_precision Az = 1;
    gradTi[0] = cos(theta)*Ar - sin(theta)*At;
    gradTi[1] = sin(theta)*Ar + cos(theta)*At;
    gradTi[2] = Az;
    */
    gitr_precision vx = p->vx[indx];
    gitr_precision vy = p->vy[indx];
    gitr_precision vz = p->vz[indx];
        vNorm = std::sqrt(vx*vx + vy*vy + vz*vz);
    p->vD[indx] = dv_ITG[2];    
	//std::cout << "gradTi Parallel " << gradTiPar << std::endl;
        //std::cout << "gradTi Parallel " << gradTi[0]<<gradTi[1]<<gradTi[2] << std::endl;
        //p->vx[indx] = p->vx[indx] +dv_ITG[0];//alpha*(gradTe[0])   
	//p->vy[indx] = p->vy[indx] +dv_ITG[1];//alpha*(gradTe[1])
	//p->vz[indx] = p->vz[indx] +dv_ITG[2];//alpha*(gradTe[2])		
        //vNorm2 = sqrt(p->vx[indx]*p->vx[indx] + p->vy[indx]*p->vy[indx] + p->vz[indx]*p->vz[indx]);
		//SFT
        gitr_precision k1 = dv_ITG[2] - dt_step*p->nu_s[indx]
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
