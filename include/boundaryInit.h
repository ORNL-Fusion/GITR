#ifndef _BOUNDARYINIT_
#define _BOUNDARYINIT_


#include "Particle.h"
#include "Boundary.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif
#include <math.h>
struct boundary_init {
    float background_Z;
    float background_amu;
    int nR_Temp;
    int nZ_Temp;
    float* TempGridr;
    float* TempGridz;
    float* ti;
    int nx;
    int nz;
    float* densityGridx;
    float* densityGridz;
    float* density;
    int nxB;
    int nzB;
    float* bfieldGridr;
    float* bfieldGridz;
    float* bfieldR;
    float* bfieldZ;
    float* bfieldT;
    
  boundary_init(float _background_Z, float _background_amu,int _nx, int _nz,float* _densityGridx, float* _densityGridz,float* _density,int _nxB,int _nzB, float* _bfieldGridr, float* _bfieldGridz,float* _bfieldR,float* _bfieldZ,float* _bfieldT,
         int _nR_Temp, int _nZ_Temp, float* _TempGridr, float* _TempGridz, float* _ti )
      : background_Z(_background_Z), background_amu(_background_amu), nx(_nx), nz(_nz), densityGridx(_densityGridx), densityGridz(_densityGridz),density(_density),nxB(_nxB),nzB(_nzB), bfieldGridr(_bfieldGridr), bfieldGridz(_bfieldGridz), bfieldR(_bfieldR), bfieldZ(_bfieldZ), bfieldT(_bfieldT),
 nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz), ti(_ti) {}

    void operator()(Boundary &b) const { 
        float midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        float midpointz = 0.5*(b.z2 - b.z1) + b.z1;
        b.density = interp2dCombined(midpointx,0.0,midpointz,nx,nz,densityGridx,densityGridz,density);
        b.ti = interp2dCombined(midpointx,0.0,midpointz,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti);
        //std::cout << "b.density " << midpointx << " " << midpointz << " " << b.density << std::endl;
        float br = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldR);   
        float bz = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldZ);
        float bt = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldT); 
        float norm_B = sqrt(br*br+bz*bz+bt*bt);
        float theta = acos((-br*b.slope_dzdx + bz)/(sqrt(br*br+bz*bz+bt*bt)*sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
 
        if (theta > 3.14159265359*0.5)
        {
            theta = acos((br*b.slope_dzdx - bz)/(sqrt(br*br+bz*bz+bt*bt)*sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
        }
        b.angle = theta*180.0/3.14159265359;
        //std::cout << "b field angle at " << midpointx << " " << midpointz << " " << b.angle << std::endl;
        b.debyeLength = sqrt(8.854187e-12*b.ti/(b.density*pow(background_Z,2)*1.60217662e-19));
        b.larmorRadius = 1.44e-4*sqrt(background_amu*b.ti/2)/(background_Z*norm_B);

    }	
};

#endif
