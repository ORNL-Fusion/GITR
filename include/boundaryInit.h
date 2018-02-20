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
    float potential;
    
    boundary_init(float _background_Z, float _background_amu,int _nx, int _nz,
          float* _densityGridx, float* _densityGridz,float* _density,int _nxB,
          int _nzB, float* _bfieldGridr, float* _bfieldGridz,float* _bfieldR,
          float* _bfieldZ,float* _bfieldT,int _nR_Temp, int _nZ_Temp,
          float* _TempGridr, float* _TempGridz, float* _ti , float _potential)

      : background_Z(_background_Z), background_amu(_background_amu), nx(_nx), nz(_nz), 
        densityGridx(_densityGridx), densityGridz(_densityGridz),density(_density),
        nxB(_nxB),nzB(_nzB), bfieldGridr(_bfieldGridr), bfieldGridz(_bfieldGridz), 
        bfieldR(_bfieldR), bfieldZ(_bfieldZ), bfieldT(_bfieldT),
        nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), 
        TempGridz(_TempGridz), ti(_ti), potential(_potential) {}

    void operator()(Boundary &b) const {
#if USE3DTETGEOM
        float midpointx = b.x1 + 0.666666667*(b.x2 + 0.5*(b.x3-b.x2)-b.x1);
        float midpointy = b.y1 + 0.666666667*(b.y2 + 0.5*(b.y3-b.y2)-b.y1);
        float midpointz = b.z1 + 0.666666667*(b.z2 + 0.5*(b.z3-b.z2)-b.z1);
#else

        float midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        float midpointy = 0.0;
        float midpointz = 0.5*(b.z2 - b.z1) + b.z1;
#endif
        b.density = interp2dCombined(midpointx,midpointy,midpointz,nx,nz,densityGridx,densityGridz,density);
        b.ti = interp2dCombined(midpointx,midpointy,midpointz,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti);
        float B[3] = {0.0,0.0,0.0};
interp2dVector(&B[0],midpointx,midpointy,midpointz,nxB,nzB,bfieldGridr,
                 bfieldGridz,bfieldR,bfieldZ,bfieldT);
        float norm_B = vectorNorm(B);
#if USE3DTETGEOM
        float surfNorm[3] = {0.0,0.0,0.0};
        b.getSurfaceNormal(surfNorm);
        float theta = acos(vectorDotProduct(B,surfNorm)/(vectorNorm(B)*vectorNorm(surfNorm)));
        if (theta > 3.14159265359*0.5)
        {
          theta = abs(theta - (3.14159265359));
        }
#else
        float theta = acos((-br*b.slope_dzdx + bz)/(sqrt(br*br+bz*bz+bt*bt)*sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
 
        if (theta > 3.14159265359*0.5)
        {
            theta = acos((br*b.slope_dzdx - bz)/(sqrt(br*br+bz*bz+bt*bt)*sqrt(b.slope_dzdx*b.slope_dzdx + 1.0)));
        }
#endif        
        b.angle = theta*180.0/3.14159265359;
        b.debyeLength = sqrt(8.854187e-12*b.ti/(b.density*pow(background_Z,2)*1.60217662e-19));
        b.larmorRadius = 1.44e-4*sqrt(background_amu*b.ti/2)/(background_Z*norm_B);
        b.flux = 0.25*b.density*sqrt(8.0*b.ti*1.602e-19/(3.1415*background_amu));
        b.impacts = 0.0;
#if BIASED_SURFACE
        b.potential = potential;
        //float cs = sqrt(2*b.ti*1.602e-19/(1.66e-27*background_amu));
        //float jsat_ion = 1.602e-19*b.density*cs;
        //b.ChildLangmuirDist = 2.0/3.0*pow(2*1.602e-19/(background_amu*1.66e-27),0.25)
        //*pow(potential,0.75)/(2.0*sqrt(3.1415*jsat_ion))*1.055e-5;
        b.ChildLangmuirDist = b.debyeLength*pow(abs(b.potential)/b.ti,0.75);
#else
        b.potential = 3.0*b.ti;    
#endif        
        //std::cout << " angle debye length and lr " <<b.angle << " " << b.debyeLength << " " << b.larmorRadius << std::endl;        
    }	
};

#endif
