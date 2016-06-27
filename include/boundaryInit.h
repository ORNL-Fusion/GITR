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
    std::vector<double>* densityGridx;
    std::vector<double>* densityGridz;
    std::vector<double>* density;
    std::vector<double>* bfieldGridr;
    std::vector<double>* bfieldGridz;
    std::vector<double>* bfieldR;
    std::vector<double>* bfieldZ;
    std::vector<double>* bfieldT;
    
  boundary_init(std::vector<double>* _densityGridx, std::vector<double>* _densityGridz,std::vector<double>* _density, std::vector<double>* _bfieldGridr, std::vector<double>* _bfieldGridz,std::vector<double>* _bfieldR,std::vector<double>* _bfieldZ,std::vector<double>* _bfieldT)
      : densityGridx(_densityGridx), densityGridz(_densityGridz),density(_density), bfieldGridr(_bfieldGridr), bfieldGridz(_bfieldGridz), bfieldR(_bfieldR), bfieldZ(_bfieldZ), bfieldT(_bfieldT) {}

    void operator()(Boundary &b) const { 
        double midpointx = 0.5*(b.x2 - b.x1)+ b.x1;
        double midpointz = 0.5*(b.z2 - b.z1) + b.z1;
        b.density = interp2d(midpointx,0.0,midpointz,densityGridx,densityGridz,density);
        std::cout << "b.density " << midpointx << " " << midpointz << " " << b.density << std::endl;
        double br = interp2d(midpointx,0.0,midpointz,bfieldGridr,bfieldGridz,bfieldR);   
        double bz = interp2d(midpointx,0.0,midpointz,bfieldGridr,bfieldGridz,bfieldZ);
        double bt = interp2d(midpointx,0.0,midpointz,bfieldGridr,bfieldGridz,bfieldT); 
        b.angle = sqrt(br*br+ bz*bz + bt*bt);
        std::cout << "b field mag at " << midpointx << " " << midpointz << " " << b.angle << std::endl;
    }	
};

#endif
