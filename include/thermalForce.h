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
    
    float dv_ITGx;
    float dv_ITGy;
    float dv_ITGz;
    float dv_ETGx;
    float dv_ETGy;
    float dv_ETGz;
    float dv_ITG[3];        
    thermalForce(Particles *_particlesPointer, float _dt, float _background_amu,
                 int _nR_gradT, int _nZ_gradT, float* _gradTGridr, float* _gradTGridz,
                 float* _gradTiR, float* _gradTiZ, float* _gradTiT, float* _gradTeR, 
                 float* _gradTeZ,float* _gradTeT,
                 int _nR_Bfield, int _nZ_Bfield,
                 float * _BfieldGridRDevicePointer,
                 float * _BfieldGridZDevicePointer,
                 float * _BfieldRDevicePointer,
                 float * _BfieldZDevicePointer,
                 float * _BfieldTDevicePointer)
        
       : particlesPointer(_particlesPointer), dt(_dt), background_amu(_background_amu),
        nR_gradT(_nR_gradT),nZ_gradT(_nZ_gradT),
        gradTGridr(_gradTGridr), gradTGridz(_gradTGridz),
        gradTiR(_gradTiR), gradTiZ(_gradTiZ),gradTiT(_gradTiT), gradTeR(_gradTeR), 
        gradTeZ(_gradTeZ),gradTeT(_gradTeT), 
        nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), 
        BfieldGridRDevicePointer(_BfieldGridRDevicePointer), 
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), 
        BfieldTDevicePointer(_BfieldTDevicePointer)
        ,dv_ITGx{0.0},dv_ITGy{0.0},dv_ITGz{0.0},dv_ETGx{0.0},dv_ETGy{0.0},dv_ETGz{0.0},
        dv_ITG{0.0,0.0,0.0}
        {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx)  
{ 
  if((particlesPointer->hitWall[indx] == 0.0) && (particlesPointer->charge[indx] > 0.0))
  {
    float MI = 1.6737236e-27;
    float kB = 1.38064852e-23;
    float KelvinPerEv = 11604.522;
    float alpha = 0.0f;
    float beta = 0.0f;
    float mu = 0.0f;
    float gradTe[3] = {0.0f};
    float gradTi[3] = {0.0f};
    //float gradTeDotB = 0.0f;
    //float gradTiDotB = 0.0f;
    float B[3] = {0.0f};
    float B_unit[3] = {0.0f};
    float gradTiPar = 0.0f;
    //float dv_ITG[3] = {0.0f};
    //float dv_ETG[3] = {0.0f};
    float x = particlesPointer->xprevious[indx];
    float y = particlesPointer->yprevious[indx];
    float z = particlesPointer->zprevious[indx];
    float Z = particlesPointer->charge[indx];
    float amu = particlesPointer->amu[indx];
    float propConst = dt/(amu*MI)*kB*KelvinPerEv;
    float fudgeFactor = 0.1;
        
    interp2dVector(&gradTi[0],x,y,z,nR_gradT,nZ_gradT,gradTGridr ,gradTGridz ,gradTiR,gradTiZ, gradTiT );
    interp2dVector(&gradTe[0],x,y,z,nR_gradT,nZ_gradT,gradTGridr ,gradTGridz ,gradTeR,gradTeZ, gradTeT );
    interp2dVector(&B[0],x,y,z,nR_Bfield,nZ_Bfield,BfieldGridRDevicePointer,
                   BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
    vectorNormalize(B,B_unit); 
    //gradTeDotB = vectorDotProduct(B_unit,gradTe); 
    //gradTiDotB = vectorDotProduct(B_unit,gradTi);
	mu = amu/(background_amu + amu);
    alpha = Z*Z*0.71;
	beta =  3.0*(mu + 7.0711*Z*Z*(1.1*powf(mu,(2.5))- 0.35*powf(mu,(1.5))) - 1.0)/(2.6 - 2*mu+ 5.4*mu*mu);
	dv_ETGx = propConst*(alpha*(gradTe[0]*B_unit[0]))*fudgeFactor;
	this->dv_ETGy = propConst*(alpha*(gradTe[1]*B_unit[1]))*fudgeFactor;
	this->dv_ETGz = propConst*(alpha*(gradTe[2]*B_unit[2]))*fudgeFactor;

	dv_ITG[0] = propConst*(beta*(gradTi[0]*B_unit[0]))*fudgeFactor;
	this->dv_ITGy = propConst*(beta*(gradTi[1]*B_unit[1]))*fudgeFactor;
	this->dv_ITGz = propConst*(beta*(gradTi[2]*B_unit[2]))*fudgeFactor;

    particlesPointer->vx[indx] = particlesPointer->vx[indx] + dv_ETGx + dv_ITGx;
	particlesPointer->vy[indx] = particlesPointer->vy[indx] + dv_ETGy + dv_ITGy;
	particlesPointer->vz[indx] = particlesPointer->vz[indx] + dv_ETGz + dv_ITGz;	
        
  }
}
};

#endif
