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

struct boundary_init {
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
    
  boundary_init(int _nx, int _nz,double* _densityGridx, double* _densityGridz,double* _density,int _nxB,int _nzB, double* _bfieldGridr, double* _bfieldGridz,double* _bfieldR,double* _bfieldZ,double* _bfieldT)
      : nx(_nx), nz(_nz), densityGridx(_densityGridx), densityGridz(_densityGridz),density(_density),nxB(_nxB),nzB(_nzB), bfieldGridr(_bfieldGridr), bfieldGridz(_bfieldGridz), bfieldR(_bfieldR), bfieldZ(_bfieldZ), bfieldT(_bfieldT) {}

    void operator()(Boundary &b) const { 
        double midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        double midpointz = 0.5*(b.z2 - b.z1) + b.z1;
        b.density = interp2dCombined(midpointx,0.0,midpointz,nx,nz,densityGridx,densityGridz,density);
        std::cout << "b.density " << midpointx << " " << midpointz << " " << b.density << std::endl;
        double br = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldR);   
        double bz = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldZ);
        double bt = interp2dCombined(midpointx,0.0,midpointz,nxB,nzB,bfieldGridr,bfieldGridz,bfieldT); 
        b.angle = sqrt(br*br+ bz*bz + bt*bt);
        std::cout << "b field mag at " << midpointx << " " << midpointz << " " << b.angle << std::endl;
    }	
};

#endif
