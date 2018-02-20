#ifndef _SURFACE_
#define _SURFACE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <cmath>
#include <math.h>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#include <stdlib.h>
#endif
CUDA_CALLABLE_MEMBER
void getBoundaryNormal(Boundary* boundaryVector,int wallIndex,float surfaceNormalVector[],float x,float y){
  #if USE3DTETGEOM > 0
           float norm_normal = boundaryVector[wallIndex].plane_norm; 
            surfaceNormalVector[0] = boundaryVector[wallIndex].a/norm_normal;
            surfaceNormalVector[1] = boundaryVector[wallIndex].b/norm_normal;
            
            surfaceNormalVector[2] = boundaryVector[wallIndex].c/norm_normal;
  #else
            float tol = 1e12;
            float norm_normal = 0.0f;
            if (boundaryVector[wallIndex].slope_dzdx == 0.0)
                {
                 surfaceNormalVector[0] = 0.0f;
                 surfaceNormalVector[1] = 0.0f;
                 surfaceNormalVector[2] = 1.0f;
                }
            else if (fabsf(boundaryVector[wallIndex].slope_dzdx)>= 0.75f*tol)
                {
                    surfaceNormalVector[0] = 1.0f;
                    surfaceNormalVector[1] = 0.0f;
                    surfaceNormalVector[2] = 0.0f;
                }
            else
                {
                    surfaceNormalVector[0] = 1.0f;
                    surfaceNormalVector[1] = 0.0f;
                    surfaceNormalVector[2] = -1.0f / (boundaryVector[wallIndex].slope_dzdx);
            norm_normal = sqrt(surfaceNormalVector[2]*surfaceNormalVector[2] + 1.0); 
            surfaceNormalVector[0] = surfaceNormalVector[0]/norm_normal;
            surfaceNormalVector[1] = surfaceNormalVector[1]/norm_normal;
            
            surfaceNormalVector[2] = surfaceNormalVector[2]/norm_normal;
                }
#if USECYLSYMM > 0 
            float theta = atan2f(y,x);
            float Sr = surfaceNormalVector[0];
            surfaceNormalVector[0] = cosf(theta)*Sr;
            surfaceNormalVector[1] = sinf(theta)*Sr;
#endif            
#endif
}
CUDA_CALLABLE_MEMBER
double screeningLength ( double Zprojectile, double Ztarget ) {
	double bohrRadius = 5.29177e-11;
	double screenLength;

	screenLength = 0.885341*bohrRadius*powf(powf(Zprojectile,(2/3)) + powf(Ztarget,(2/3)),(-1/2));

	return screenLength;
}

CUDA_CALLABLE_MEMBER
double stoppingPower (Particles * particles,int indx, double Mtarget, double Ztarget, double screenLength) {
	        double E0;
            double Q = 1.60217662e-19;
		double ke2 = 14.4e-10;
		double reducedEnergy;
	double stoppingPower;

	E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/Q;
	reducedEnergy = E0*(Mtarget/(particles->amu[indx]+Mtarget))*(screenLength/(particles->Z[indx]*Ztarget*ke2));
	stoppingPower = 0.5*log(1.0 + 1.2288*reducedEnergy)/(reducedEnergy + 0.1728*sqrt(reducedEnergy) + 0.008*pow(reducedEnergy, 0.1504));

	return stoppingPower;	
}

struct erosion { 
    Particles *particles;
    const double dt;

    erosion(Particles *_particles, double _dt) : particles(_particles), dt(_dt) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 
	double screenLength;
	double stopPower;
	double q = 18.6006;
	double lambda = 2.2697;
    double mu = 3.1273;
	double Eth = 24.9885;
	double Y0;
	double Ztarget = 74.0;
	double Mtarget = 183.84;
	double term;
	double E0;

	screenLength = screeningLength(particles->Z[indx], Ztarget);
	stopPower = stoppingPower(particles,indx, Mtarget, Ztarget, screenLength); 
	E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/1.60217662e-19;
	term = pow((E0/Eth - 1),mu);
	Y0 = q*stopPower*term/(lambda + term);
    	}
     
};

struct reflection {
    Particles * particles;
    const double dt;
    int nLines;
    Boundary * boundaryVector;
    int nE_sputtRefCoeff;
    int nA_sputtRefCoeff;
    float* A_sputtRefCoeff;
    float* Elog_sputtRefCoeff;
    float* spyl_surfaceModel;
    float* rfyl_surfaceModel;
    int nE_sputtRefDistOut; 
    int nA_sputtRefDistOut;
    int nE_sputtRefDistIn;
    int nA_sputtRefDistIn;
    float* E_sputtRefDistIn;
    float* A_sputtRefDistIn;
    float* E_sputtRefDistOut;
    float* A_sputtRefDistOut;
    float* energyDistGrid01;
    float* angleDistGrid01;
    float* EDist_CDF_Y_regrid;
    float* ADist_CDF_Y_regrid;
    float* EDist_CDF_R_regrid;
    float* ADist_CDF_R_regrid;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    reflection(Particles* _particles, double _dt,
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nLines,Boundary * _boundaryVector,
    int _nE_sputtRefCoeff,
    int _nA_sputtRefCoeff,
    float* _A_sputtRefCoeff,
    float* _Elog_sputtRefCoeff,
    float* _spyl_surfaceModel,
    float* _rfyl_surfaceModel,
    int _nE_sputtRefDistOut,
    int _nA_sputtRefDistOut,
    int _nE_sputtRefDistIn,
    int _nA_sputtRefDistIn,
    float* _E_sputtRefDistIn,
    float* _A_sputtRefDistIn,
    float* _E_sputtRefDistOut,
    float* _A_sputtRefDistOut,
    float* _energyDistGrid01,
    float* _angleDistGrid01,
    float* _EDist_CDF_Y_regrid,
    float* _ADist_CDF_Y_regrid, 
    float* _EDist_CDF_R_regrid,
    float* _ADist_CDF_R_regrid) :
    particles(_particles), dt(_dt), state(_state),nLines(_nLines),boundaryVector(_boundaryVector),
    nE_sputtRefCoeff(_nE_sputtRefCoeff),nA_sputtRefCoeff(_nA_sputtRefCoeff),
    A_sputtRefCoeff(_A_sputtRefCoeff),
    Elog_sputtRefCoeff(_Elog_sputtRefCoeff),
    spyl_surfaceModel(_spyl_surfaceModel),
    rfyl_surfaceModel(_rfyl_surfaceModel),
    nE_sputtRefDistOut(_nE_sputtRefDistOut),
    nA_sputtRefDistOut(_nA_sputtRefDistOut),
    nE_sputtRefDistIn(_nE_sputtRefDistIn),
    nA_sputtRefDistIn(_nA_sputtRefDistIn),
    E_sputtRefDistIn(_E_sputtRefDistIn),
    A_sputtRefDistIn(_A_sputtRefDistIn),
    E_sputtRefDistOut(_E_sputtRefDistOut),
    A_sputtRefDistOut(_A_sputtRefDistOut),
    energyDistGrid01(_energyDistGrid01),
    angleDistGrid01(_angleDistGrid01),
    EDist_CDF_Y_regrid(_EDist_CDF_Y_regrid),
    ADist_CDF_Y_regrid(_ADist_CDF_Y_regrid),
    EDist_CDF_R_regrid(_EDist_CDF_R_regrid),
    ADist_CDF_R_regrid(_ADist_CDF_R_regrid){}

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const {
    
  if(particles->hitWall[indx] == 1.0)
  {   
    float E0 = 0.0;
    float thetaImpact = 0.0;
    float particleTrackVector[3] = {0.0f};
    float surfaceNormalVector[3] = {0.0f};
    float vSampled[3] = {0.0f};
    float norm_part = 0.0;
    int signPartDotNormal=0; 
    float partDotNormal = 0.0;
    float Enew = 0.0f;
    float angleSample = 0.0f;
    int wallIndex = 0;
    float tol = 1e12;
    float Sr = 0.0;
    float St = 0.0;
    float Y0 = 0.0;
    float R0 = 0.0;
    float totalYR=0.0;
    float newWeight=0.0;
    int wallHit = particles->wallHit[indx];
    float eInterpVal=0.0;
    float aInterpVal=0.0;

    E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/1.60217662e-19;
    particleTrackVector[0] = particles->vx[indx];
    particleTrackVector[1] = particles->vy[indx];
    particleTrackVector[2] = particles->vz[indx];
    wallIndex = particles->wallIndex[indx];
    boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector);
    //#if USECYLSYMM > 0 
    //  float theta = atan2f(particles->yprevious[indx],particles->xprevious[indx]);
    //  Sr = surfaceNormalVector[0];
    //  surfaceNormalVector[0] = cosf(theta)*Sr;
    //  surfaceNormalVector[1] = sinf(theta)*Sr;
    //#endif            
            //std::cout << "velocities " << particles->vx[indx] << " " << particles->vy[indx] << " " << particles->vz[indx] << std::endl;
            //std::cout << "surface norm " << surfaceNormalVector[0] << " " << surfaceNormalVector[1] << " " << surfaceNormalVector[2] << std::endl;
    norm_part = sqrt(particleTrackVector[0]*particleTrackVector[0] + particleTrackVector[1]*particleTrackVector[1] + particleTrackVector[2]*particleTrackVector[2]);
            //          norm_normal = sqrt(surfaceNormalVectorRotated[0]*surfaceNormalVectorRotated[0] + surfaceNormalVectorRotated[1]*surfaceNormalVectorRotated[1] + surfaceNormalVectorRotated[2]*surfaceNormalVectorRotated[2]);
    //        std::cout << "norm of particle track " << norm_part << std::endl;
    //        std::cout << "norm of surface normal " << norm_normal << std::endl;
            particleTrackVector[0] = particleTrackVector[0]/norm_part;
            particleTrackVector[1] = particleTrackVector[1]/norm_part;
            particleTrackVector[2] = particleTrackVector[2]/norm_part;

            partDotNormal = vectorDotProduct(particleTrackVector,surfaceNormalVector);
            thetaImpact = acos(partDotNormal);
            if (thetaImpact > 3.14159265359*0.5)
            {
              thetaImpact = abs(thetaImpact - (3.14159265359));
            }
            thetaImpact = thetaImpact*180.0/3.14159265359;
            signPartDotNormal = sgn(partDotNormal);
            if( boundaryVector[wallHit].Z > 0.0)
            {
                Y0 = interp2d(thetaImpact,log10(E0),nA_sputtRefCoeff,
                    nE_sputtRefCoeff,A_sputtRefCoeff,
                    Elog_sputtRefCoeff,spyl_surfaceModel);
                R0 = interp2d(thetaImpact,log10(E0),nA_sputtRefCoeff, 
                    nE_sputtRefCoeff,A_sputtRefCoeff,
                    Elog_sputtRefCoeff,rfyl_surfaceModel);
            }
            else
            {
                Y0 = 0.0;
                R0 = 0.0;
            }
            totalYR=Y0+R0;

            //std::cout << "Energy angle yield " << E0 << " " << thetaImpact << " " << Y0 << std::endl; 
            #if PARTICLESEEDS > 0
              #ifdef __CUDACC__
                float r7 = curand_uniform(&state[indx]);
                float r8 = curand_uniform(&state[indx]);
                float r9 = curand_uniform(&state[indx]);
                float r10 = curand_uniform(&state[indx]);
              #else
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                float r7 = dist(state[indx]);
                float r8 = dist(state[indx]);
                float r9 = dist(state[indx]);
                float r10 = dist(state[indx]);
              #endif

            #else
              #if __CUDACC__
                float r7 = curand_uniform(&state[6]);
                float r8 = curand_uniform(&state[7]);
                float r9 = curand_uniform(&state[8]);
              #else
                std::uniform_real_distribution<float> dist(0.0, 1.0);
                float r7=dist(state[6]);
                float r8=dist(state[7]);
                float r9=dist(state[8]);
              #endif
                //float r7 = 0.0;
            #endif
            //particle either reflects or deposits
            float sputtProb = Y0/totalYR;
            if(r7 > sputtProb) //reflects
            {
                  aInterpVal = interp3d (r8,thetaImpact,log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                    angleDistGrid01,A_sputtRefDistIn,
                                    E_sputtRefDistIn,ADist_CDF_R_regrid);
                   eInterpVal = interp3d ( r9,thetaImpact,log10(E0),
                           nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                         energyDistGrid01,A_sputtRefDistIn,
                                         E_sputtRefDistIn,EDist_CDF_R_regrid );
                   newWeight=R0/(1.0f-sputtProb);
            }
            else //sputters
            {
                  aInterpVal = interp3d(r8,thetaImpact,log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                          angleDistGrid01,A_sputtRefDistIn,
                          E_sputtRefDistIn,ADist_CDF_Y_regrid);
                  eInterpVal = interp3d(r9,thetaImpact,log10(E0),
                           nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                           energyDistGrid01,A_sputtRefDistIn,
                           E_sputtRefDistIn,EDist_CDF_Y_regrid);
                   newWeight=Y0/sputtProb;
            }
                //deposit on surface
            //#if USESURFACEMODEL == 1
            //  #if USE_CUDA > 0
            //    atomicAdd(&boundaryVector[wallHit].impacts, particles->weight[indx]);
            //  #else
            //    boundaryVector[wallHit].impacts = boundaryVector[wallHit].impacts +  particles->weight[indx];
            //  #endif
            //#endif
            //float deposited = 0.0; 
            //   if(Y0 <= 1.0)
            //   {
            //     deposited = particles->weight[indx]*(1.0-Y0);
            //   }
            //#if USE_CUDA > 0
            //    atomicAdd(&boundaryVector[wallHit].redeposit, deposited);
            //#else
            //   boundaryVector[wallHit].redeposit = boundaryVector[wallHit].redeposit +deposited;
            //#endif
                //reflect with weight and new initial conditions
                if( boundaryVector[wallHit].Z > 0.0)
                {
                particles->weight[indx] = particles->weight[indx]*newWeight;
                particles->hitWall[indx] = 0.0;
                particles->charge[indx] = 0.0;
                float V0 = sqrt(2*eInterpVal*1.602e-19/(particles->amu[indx]*1.66e-27));
    vSampled[0] = V0*sin(aInterpVal*3.1415/180)*cos(2.0*3.1415*r10);
    vSampled[1] = V0*sin(aInterpVal*3.1415/180)*sin(2.0*3.1415*r10);
    vSampled[2] = V0*cos(aInterpVal*3.1415/180);
    boundaryVector[wallHit].transformToSurface(vSampled);
    particles->vx[indx] = -signPartDotNormal*vSampled[0];
    particles->vy[indx] = -signPartDotNormal*vSampled[1];
    particles->vz[indx] = -signPartDotNormal*vSampled[2];
    particles->test[indx] = -signPartDotNormal*vSampled[0];
    particles->test0[indx] = -signPartDotNormal*vSampled[1];
    particles->test1[indx] = -signPartDotNormal*vSampled[2];

    particles->xprevious[indx] = particles->x[indx] + -signPartDotNormal*surfaceNormalVector[0]*1e-6;
    particles->yprevious[indx] = particles->y[indx] + -signPartDotNormal*surfaceNormalVector[1]*1e-6;
    particles->zprevious[indx] = particles->z[indx] + -signPartDotNormal*surfaceNormalVector[2]*1e-6;
            }
                else
                {
                    particles->hitWall[indx] = 2.0;
                }
                  }
}
};
#endif
