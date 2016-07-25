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
    double background_Z;
    double background_amu;
    int nR_Temp;
    int nZ_Temp;
    double* TempGridr;
    double* TempGridz;
    double* ti;
    int nx;
    int nz;
    double* densityGridx;
    double* densityGridz;
    double* density;
    int nxB;
    int nzB;
    double* bfieldGridr;
    double* bfieldGridz;
    double* bfieldR;
    double* bfieldZ;
    double* bfieldT;
    
  boundary_init(double _background_Z, double _background_amu,int _nx, int _nz,double* _densityGridx, double* _densityGridz,double* _density,int _nxB,int _nzB, double* _bfieldGridr, double* _bfieldGridz,double* _bfieldR,double* _bfieldZ,double* _bfieldT,
         int _nR_Temp, int _nZ_Temp, double* _TempGridr, double* _TempGridz, double* _ti )
      : background_Z(_background_Z), background_amu(_background_amu), nx(_nx), nz(_nz), densityGridx(_densityGridx), densityGridz(_densityGridz),density(_density),nxB(_nxB),nzB(_nzB), bfieldGridr(_bfieldGridr), bfieldGridz(_bfieldGridz), bfieldR(_bfieldR), bfieldZ(_bfieldZ), bfieldT(_bfieldT),
 nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz), ti(_ti) {}

    void operator()(Boundary &b) const { 
        double midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        double midpointz = 0.5*(b.z2 - b.z1) + b.z1;
        b.density = interp2dCombined(midpointx,0.0,midpointz,nx,nz,densityGridx,densityGridz,density);
        b.ti = interp2dCombined(midpointx,0.0,midpointz,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti);
//        std::cout << "b.density " << midpointx << " " << midpointz << " " << b.density << std::endl;
        double br = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldR);   
        double bz = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldZ);
        double bt = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldT); 
        double norm_B = sqrt(br*br+bz*bz+bt*bt);
        double theta = acos((-br*b.slope_dzdx + bz)/(sqrt(br*br+bz*bz+bt*bt)*sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
 
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
