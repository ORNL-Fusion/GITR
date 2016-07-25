#ifndef _IONIZE_
#define _IONIZE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particle.h"
#ifdef __CUDACC__
#include <thrust/random.h>
#endif

#ifdef __GNUC__ 
#include <random>
#include <stdlib.h>
#endif

#include "interpRateCoeff.hpp"

struct ionize { 
    int nR_Dens;
    int nZ_Dens;
    double* DensGridr;
    double* DensGridz;
    double* ne;
    int nR_Temp;
    int nZ_Temp;
    double* TempGridr;
    double* TempGridz;
    double* te;
    int nTemperaturesIonize;
    int nDensitiesIonize;
    double* gridDensity_Ionization;
    double* gridTemperature_Ionization;
    double* rateCoeff_Ionization;
    const double dt;
        ionize(double _dt,int _nR_Dens,int _nZ_Dens,double* _DensGridr,
    double* _DensGridz,double* _ne,int _nR_Temp, int _nZ_Temp,
    double* _TempGridr, double* _TempGridz,double* _te,int _nTemperaturesIonize,
    int _nDensitiesIonize,double* _gridTemperature_Ionization,double* _gridDensity_Ionization,
    double* _rateCoeff_Ionization
              ) : dt(_dt), nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),ne(_ne),
        nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz), te(_te),
   nTemperaturesIonize(_nTemperaturesIonize), nDensitiesIonize(_nDensitiesIonize), gridTemperature_Ionization(_gridTemperature_Ionization), gridDensity_Ionization(_gridDensity_Ionization), rateCoeff_Ionization(_rateCoeff_Ionization) {}
    
        CUDA_CALLABLE_MEMBER_DEVICE 
                void operator()(Particle &p) const { 
	if(p.hitWall == 0.0){        
	double tion = interpRateCoeff2d ( p.charge, p.x, p.y, p.z,nR_Temp,nZ_Temp, TempGridr,TempGridz,te,DensGridr,DensGridz, ne,nTemperaturesIonize,nDensitiesIonize,gridTemperature_Ionization,gridDensity_Ionization,rateCoeff_Ionization );	
    double P1 = 1-exp(-dt/tion);
	
	#ifdef __CUDACC__
	double r1 = curand_uniform(&p.streams[0]);
	#else
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	double r1=dist(p.streams[0]);
	#endif

	if(r1 <= P1)
	{
		p.charge = p.charge+1;
	} 
	}	

	} 
};

#endif
