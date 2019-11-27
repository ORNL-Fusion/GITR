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
#include "Surfaces.h"
#include <cmath>
#include "boris.h"

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
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
            else if (std::abs(boundaryVector[wallIndex].slope_dzdx)>= 0.75f*tol)
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
            norm_normal = std::sqrt(surfaceNormalVector[2]*surfaceNormalVector[2] + 1.0); 
            surfaceNormalVector[0] = surfaceNormalVector[0]/norm_normal;
            surfaceNormalVector[1] = surfaceNormalVector[1]/norm_normal;
            
            surfaceNormalVector[2] = surfaceNormalVector[2]/norm_normal;
                }
#if USECYLSYMM > 0 
            float theta = std::atan2(y,x);
            float Sr = surfaceNormalVector[0];
            surfaceNormalVector[0] = std::cos(theta)*Sr;
            surfaceNormalVector[1] = std::sin(theta)*Sr;
#endif            
#endif
}
CUDA_CALLABLE_MEMBER
double screeningLength ( double Zprojectile, double Ztarget ) {
	double bohrRadius = 5.29177e-11;
	double screenLength;

	screenLength = 0.885341*bohrRadius*std::pow(std::pow(Zprojectile,(2/3)) + std::pow(Ztarget,(2/3)),(-1/2));

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
	stoppingPower = 0.5*std::log(1.0 + 1.2288*reducedEnergy)/(reducedEnergy + 0.1728*std::sqrt(reducedEnergy) + 0.008*std::pow(reducedEnergy, 0.1504));

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
	term = std::pow((E0/Eth - 1),mu);
	Y0 = q*stopPower*term/(lambda + term);
    	}
     
};

struct reflection {
    Particles * particles;
    const double dt;
    int nLines;
    Boundary * boundaryVector;
    Surfaces * surfaces;
    int nE_sputtRefCoeff;
    int nA_sputtRefCoeff;
    float* A_sputtRefCoeff;
    float* Elog_sputtRefCoeff;
    float* spyl_surfaceModel;
    float* rfyl_surfaceModel;
    int nE_sputtRefDistOut; 
    int nE_sputtRefDistOutRef; 
    int nA_sputtRefDistOut;
    int nE_sputtRefDistIn;
    int nA_sputtRefDistIn;
    float* E_sputtRefDistIn;
    float* A_sputtRefDistIn;
    float* E_sputtRefDistOut;
    float* E_sputtRefDistOutRef;
    float* A_sputtRefDistOut;
    float* energyDistGrid01;
    float* energyDistGrid01Ref;
    float* angleDistGrid01;
    float* EDist_CDF_Y_regrid;
    float* ADist_CDF_Y_regrid;
    float* EDist_CDF_R_regrid;
    float* ADist_CDF_R_regrid;
    int nEdist;
    float E0dist;
    float Edist;
    int nAdist;
    float A0dist;
    float Adist;
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
            Surfaces * _surfaces,
    int _nE_sputtRefCoeff,
    int _nA_sputtRefCoeff,
    float* _A_sputtRefCoeff,
    float* _Elog_sputtRefCoeff,
    float* _spyl_surfaceModel,
    float* _rfyl_surfaceModel,
    int _nE_sputtRefDistOut,
    int _nE_sputtRefDistOutRef,
    int _nA_sputtRefDistOut,
    int _nE_sputtRefDistIn,
    int _nA_sputtRefDistIn,
    float* _E_sputtRefDistIn,
    float* _A_sputtRefDistIn,
    float* _E_sputtRefDistOut,
    float* _E_sputtRefDistOutRef,
    float* _A_sputtRefDistOut,
    float* _energyDistGrid01,
    float* _energyDistGrid01Ref,
    float* _angleDistGrid01,
    float* _EDist_CDF_Y_regrid,
    float* _ADist_CDF_Y_regrid, 
    float* _EDist_CDF_R_regrid,
    float* _ADist_CDF_R_regrid,
    int _nEdist,
    float _E0dist,
    float _Edist,
    int _nAdist,
    float _A0dist,
    float _Adist) :
particles(_particles),
                             dt(_dt),
                             nLines(_nLines),
                             boundaryVector(_boundaryVector),
                             surfaces(_surfaces),
                             nE_sputtRefCoeff(_nE_sputtRefCoeff),
                             nA_sputtRefCoeff(_nA_sputtRefCoeff),
                             A_sputtRefCoeff(_A_sputtRefCoeff),
                             Elog_sputtRefCoeff(_Elog_sputtRefCoeff),
                             spyl_surfaceModel(_spyl_surfaceModel),
                             rfyl_surfaceModel(_rfyl_surfaceModel),
                             nE_sputtRefDistOut(_nE_sputtRefDistOut),
                             nE_sputtRefDistOutRef(_nE_sputtRefDistOutRef),
                             nA_sputtRefDistOut(_nA_sputtRefDistOut),
                             nE_sputtRefDistIn(_nE_sputtRefDistIn),
                             nA_sputtRefDistIn(_nA_sputtRefDistIn),
                             E_sputtRefDistIn(_E_sputtRefDistIn),
                             A_sputtRefDistIn(_A_sputtRefDistIn),
                             E_sputtRefDistOut(_E_sputtRefDistOut),
                             E_sputtRefDistOutRef(_E_sputtRefDistOutRef),
                             A_sputtRefDistOut(_A_sputtRefDistOut),
                             energyDistGrid01(_energyDistGrid01),
                             energyDistGrid01Ref(_energyDistGrid01Ref),
                             angleDistGrid01(_angleDistGrid01),
                             EDist_CDF_Y_regrid(_EDist_CDF_Y_regrid),
                             ADist_CDF_Y_regrid(_ADist_CDF_Y_regrid),
                             EDist_CDF_R_regrid(_EDist_CDF_R_regrid),
                             ADist_CDF_R_regrid(_ADist_CDF_R_regrid),
                             nEdist(_nEdist),
                             E0dist(_E0dist),
                             Edist(_Edist),
                             nAdist(_nAdist),
                             A0dist(_A0dist),
                             Adist(_Adist),
                             state(_state) {
  }

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const {
    
    if (particles->hitWall[indx] == 1.0) {
      float E0 = 0.0;
      float thetaImpact = 0.0;
      float particleTrackVector[3] = {0.0f};
      float surfaceNormalVector[3] = {0.0f};
      float vSampled[3] = {0.0f};
      float norm_part = 0.0;
      int signPartDotNormal = 0;
      float partDotNormal = 0.0;
      float Enew = 0.0f;
      float angleSample = 0.0f;
      int wallIndex = 0;
      float tol = 1e12;
      float Sr = 0.0;
      float St = 0.0;
      float Y0 = 0.0;
      float R0 = 0.0;
      float totalYR = 0.0;
      float newWeight = 0.0;
      int wallHit = particles->wallHit[indx];
      int surfaceHit = boundaryVector[wallHit].surfaceNumber;
      int surface = boundaryVector[wallHit].surface;
      if (wallHit > 260)
        wallHit = 260;
      if (wallHit < 0)
        wallHit = 0;
      if (surfaceHit > 260)
        surfaceHit = 260;
      if (surfaceHit < 0)
        surfaceHit = 0;
      if (surface > 260)
        surface = 260;
      if (surface < 0)
        surface = 0;
      float eInterpVal = 0.0;
      float aInterpVal = 0.0;
      float weight = particles->weight[indx];
      float vx = particles->vx[indx];
      float vy = particles->vy[indx];
      float vz = particles->vz[indx];
#if FLUX_EA > 0
      float dEdist = (Edist - E0dist) / static_cast<float>(nEdist);
      float dAdist = (Adist - A0dist) / static_cast<float>(nAdist);
      int AdistInd = 0;
      int EdistInd = 0;
#endif
      particles->firstCollision[indx] = 1;
      particleTrackVector[0] = vx;
      particleTrackVector[1] = vy;
      particleTrackVector[2] = vz;
      norm_part = std::sqrt(particleTrackVector[0] * particleTrackVector[0] + particleTrackVector[1] * particleTrackVector[1] + particleTrackVector[2] * particleTrackVector[2]);
      E0 = 0.5 * particles->amu[indx] * 1.6737236e-27 * (norm_part * norm_part) / 1.60217662e-19;
      if (E0 > 1000.0)
        E0 = 990.0;
      //std::cout << "Particle hit wall with energy " << E0 << std::endl;
      //std::cout << "Particle hit wall with v " << vx << " " << vy << " " << vz<< std::endl;
      //std::cout << "Particle amu norm_part " << particles->amu[indx] << " " << vy << " " << vz<< std::endl;
      wallIndex = particles->wallIndex[indx];
      boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector, particles->y[indx], particles->x[indx]);
      particleTrackVector[0] = particleTrackVector[0] / norm_part;
      particleTrackVector[1] = particleTrackVector[1] / norm_part;
      particleTrackVector[2] = particleTrackVector[2] / norm_part;

      partDotNormal = vectorDotProduct(particleTrackVector, surfaceNormalVector);
      thetaImpact = std::acos(partDotNormal);
      if (thetaImpact > 3.14159265359 * 0.5) {
        thetaImpact = std::abs(thetaImpact - (3.14159265359));
      }
      thetaImpact = thetaImpact * 180.0 / 3.14159265359;
      if (thetaImpact < 0.0)
        thetaImpact = 0.0;
      signPartDotNormal = std::copysign(1.0,partDotNormal);
      if (E0 == 0.0) {
        thetaImpact = 0.0;
      }
      if (boundaryVector[wallHit].Z > 0.0) {
        Y0 = interp2d(thetaImpact, std::log10(E0), nA_sputtRefCoeff,
                      nE_sputtRefCoeff, A_sputtRefCoeff,
                      Elog_sputtRefCoeff, spyl_surfaceModel);
        R0 = interp2d(thetaImpact, std::log10(E0), nA_sputtRefCoeff,
                      nE_sputtRefCoeff, A_sputtRefCoeff,
                      Elog_sputtRefCoeff, rfyl_surfaceModel);
      } else {
        Y0 = 0.0;
        R0 = 0.0;
      }
      //std::cout << "Particle " << indx << " struck surface with energy and angle " << E0 << " " << thetaImpact << std::endl;
      //std::cout << " resulting in Y0 and R0 of " << Y0 << " " << R0 << std::endl;
      totalYR = Y0 + R0;
//if(particles->test[indx] == 0.0)
//{
//    particles->test[indx] = 1.0;
//    particles->test0[indx] = E0;
//    particles->test1[indx] = thetaImpact;
//    particles->test2[indx] = Y0;
//    particles->test3[indx] = R0;
//}
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
	    int didReflect = 0;
            if(totalYR > 0.0)
            {
            if(r7 > sputtProb) //reflects
            {
	          didReflect = 1;
                  aInterpVal = interp3d (r8,thetaImpact,std::log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                    angleDistGrid01,A_sputtRefDistIn,
                                    E_sputtRefDistIn,ADist_CDF_R_regrid);
                   eInterpVal = interp3d ( r9,thetaImpact,std::log10(E0),
                           nE_sputtRefDistOutRef,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                         energyDistGrid01Ref,A_sputtRefDistIn,
                                         E_sputtRefDistIn,EDist_CDF_R_regrid );
                   //newWeight=(R0/(1.0f-sputtProb))*weight;
		   newWeight = weight*(totalYR);
    #if FLUX_EA > 0
              EdistInd = std::floor((eInterpVal-E0dist)/dEdist);
              AdistInd = std::floor((aInterpVal-A0dist)/dAdist);
              if((EdistInd >= 0) && (EdistInd < nEdist) && 
                 (AdistInd >= 0) && (AdistInd < nAdist))
              {
            #if USE_CUDA > 0
                  atomicAdd(&surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd],newWeight);
            #else      
                  surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                    surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
            #endif
               }
	       #endif
                  if(surface > 0)
                {

            #if USE_CUDA > 0
                    atomicAdd(&surfaces->grossDeposition[surfaceHit],weight*(1.0-R0));
            #else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit]+weight*(1.0-R0);
            #endif
                }
            }
            else //sputters
            {
                  aInterpVal = interp3d(r8,thetaImpact,std::log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                          angleDistGrid01,A_sputtRefDistIn,
                          E_sputtRefDistIn,ADist_CDF_Y_regrid);
                  eInterpVal = interp3d(r9,thetaImpact,std::log10(E0),
                           nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                           energyDistGrid01,A_sputtRefDistIn,E_sputtRefDistIn,EDist_CDF_Y_regrid);
            //if(particles->test[indx] == 0.0)
            //{
            //    particles->test[indx] = 1.0;
            //    particles->test0[indx] = aInterpVal;
            //    particles->test1[indx] = eInterpVal;
            //    particles->test2[indx] = r8;
            //    particles->test3[indx] = r9;
            //}            
                //std::cout << " particle sputters with " << eInterpVal << aInterpVal <<  std::endl;
                  //newWeight=(Y0/sputtProb)*weight;
		  newWeight=weight*totalYR;
    #if FLUX_EA > 0
              EdistInd = std::floor((eInterpVal-E0dist)/dEdist);
              AdistInd = std::floor((aInterpVal-A0dist)/dAdist);
              if((EdistInd >= 0) && (EdistInd < nEdist) && 
                 (AdistInd >= 0) && (AdistInd < nAdist))
              {
                //std::cout << " particle sputters with " << EdistInd << AdistInd <<  std::endl;
            #if USE_CUDA > 0
                  atomicAdd(&surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd],newWeight);
            #else      
                  surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                    surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
              #endif 
              }
	       #endif
                  if(sputtProb == 0.0) newWeight = 0.0;
                   //std::cout << " particle sputtered with newWeight " << newWeight << std::endl;
                  if(surface > 0)
                {

            #if USE_CUDA > 0
                    atomicAdd(&surfaces->grossDeposition[surfaceHit],weight*(1.0-R0));
                    atomicAdd(&surfaces->grossErosion[surfaceHit],newWeight);
                    atomicAdd(&surfaces->aveSputtYld[surfaceHit],Y0);
                    if(weight > 0.0)
                    {
                        atomicAdd(&surfaces->sputtYldCount[surfaceHit],1);
                    }
            #else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit]+weight*(1.0-R0);
                    surfaces->grossErosion[surfaceHit] = surfaces->grossErosion[surfaceHit] + newWeight;
                    surfaces->aveSputtYld[surfaceHit] = surfaces->aveSputtYld[surfaceHit] + Y0;
                    surfaces->sputtYldCount[surfaceHit] = surfaces->sputtYldCount[surfaceHit] + 1;
            #endif
                }
            }
            //std::cout << "eInterpValYR " << eInterpVal << std::endl; 
            }
            else
            {       newWeight = 0.0;
                    particles->hitWall[indx] = 2.0;
                  if(surface > 0)
                {
            #if USE_CUDA > 0
                    atomicAdd(&surfaces->grossDeposition[surfaceHit],weight);
            #else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit]+weight;
            #endif
	        }
            //std::cout << "eInterpValYR_not " << eInterpVal << std::endl; 
            }
            //std::cout << "eInterpVal " << eInterpVal << std::endl; 
	    if(eInterpVal <= 0.0)
            {       newWeight = 0.0;
                    particles->hitWall[indx] = 2.0;
                  if(surface > 0)
                {
		    if(didReflect)
		    {
            #if USE_CUDA > 0
                    atomicAdd(&surfaces->grossDeposition[surfaceHit],weight);
            #else
                    surfaces->grossDeposition[surfaceHit] = surfaces->grossDeposition[surfaceHit]+weight;
            #endif
	            }
		}
            }
            //if(particles->test[indx] == 1.0)
            //{
            //    particles->test3[indx] = eInterpVal;
            //    particles->test[indx] = 2.0;
            //}
                //deposit on surface
            if(surface)
            {
            #if USE_CUDA > 0
                atomicAdd(&surfaces->sumWeightStrike[surfaceHit],weight);
                atomicAdd(&surfaces->sumParticlesStrike[surfaceHit],1);
            #else
                surfaces->sumWeightStrike[surfaceHit] =surfaces->sumWeightStrike[surfaceHit] +weight;
                surfaces->sumParticlesStrike[surfaceHit] = surfaces->sumParticlesStrike[surfaceHit]+1;
              //boundaryVector[wallHit].impacts = boundaryVector[wallHit].impacts +  particles->weight[indx];
            #endif
            #if FLUX_EA > 0
                EdistInd = std::floor((E0-E0dist)/dEdist);
                AdistInd = std::floor((thetaImpact-A0dist)/dAdist);
              
	        if((EdistInd >= 0) && (EdistInd < nEdist) && 
                  (AdistInd >= 0) && (AdistInd < nAdist))
                {
#if USE_CUDA > 0
                    atomicAdd(&surfaces->energyDistribution[surfaceHit*nEdist*nAdist + 
                                               EdistInd*nAdist + AdistInd], weight);
#else

                    surfaces->energyDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                    surfaces->energyDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  weight;
#endif
        }
#endif
      }
      //reflect with weight and new initial conditions
      //std::cout << "particle wall hit Z and nwweight " << boundaryVector[wallHit].Z << " " << newWeight << std::endl;
      if (boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0)
      //if(newWeight > 0.0)
      {
        particles->weight[indx] = newWeight;
        particles->hitWall[indx] = 0.0;
        particles->charge[indx] = 0.0;
        float V0 = std::sqrt(2 * eInterpVal * 1.602e-19 / (particles->amu[indx] * 1.66e-27));
        particles->newVelocity[indx] = V0;
        vSampled[0] = V0 * std::sin(aInterpVal * 3.1415 / 180) * std::cos(2.0 * 3.1415 * r10);
        vSampled[1] = V0 * std::sin(aInterpVal * 3.1415 / 180) * std::sin(2.0 * 3.1415 * r10);
        vSampled[2] = V0 * std::cos(aInterpVal * 3.1415 / 180);
        boundaryVector[wallHit].transformToSurface(vSampled, particles->y[indx], particles->x[indx]);
        //float rr = std::sqrt(particles->x[indx]*particles->x[indx] + particles->y[indx]*particles->y[indx]);
        //if (particles->z[indx] < -4.1 && -signPartDotNormal*vSampled[0] > 0.0)
        //{
        //  std::cout << "particle index " << indx  << std::endl;
        //  std::cout << "aInterpVal" << aInterpVal  << std::endl;
        //  std::cout << "Surface Normal" << surfaceNormalVector[0] << " " << surfaceNormalVector[1] << " " << surfaceNormalVector[2] << std::endl;
        //  std::cout << "signPartDotNormal " << signPartDotNormal << std::endl;
        //  std::cout << "Particle hit wall with v " << vx << " " << vy << " " << vz<< std::endl;
        //  std::cout << "vSampled " << vSampled[0] << " " << vSampled[1] << " " << vSampled[2] << std::endl;
        //  std::cout << "Final transform" << -signPartDotNormal*vSampled[0] << " " << -signPartDotNormal*vSampled[1] << " " << -signPartDotNormal*vSampled[2] << std::endl;
        //  std::cout << "Position of particle0 " << particles->xprevious[indx] << " " << particles->yprevious[indx] << " " << particles->zprevious[indx] << std::endl;
        //  std::cout << "Position of particle " << particles->x[indx] << " " << particles->y[indx] << " " << particles->z[indx] << std::endl;
        //  }
        particles->vx[indx] = -static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * vSampled[0];
        particles->vy[indx] = -static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * vSampled[1];
        particles->vz[indx] = -static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * vSampled[2];
        //        //if(particles->test[indx] == 0.0)
        //        //{
        //        //    particles->test[indx] = 1.0;
        //        //    particles->test0[indx] = aInterpVal;
        //        //    particles->test1[indx] = eInterpVal;
        //        //    particles->test2[indx] = V0;
        //        //    particles->test3[indx] = vSampled[2];
        //        //}

        particles->xprevious[indx] = particles->x[indx] - static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[0] * 1e-4;
        particles->yprevious[indx] = particles->y[indx] - static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[1] * 1e-4;
      particles->zprevious[indx] = particles->z[indx] - static_cast<float>(boundaryVector[wallHit].inDir) * surfaceNormalVector[2] * 1e-4;
        //std::cout << "New vel " << particles->vx[indx] << " " << particles->vy[indx] << " " << particles->vz[indx] << std::endl;
        //std::cout << "New pos " << particles->xprevious[indx] << " " << particles->yprevious[indx] << " " << particles->zprevious[indx] << std::endl;
        //if(particles->test[indx] == 0.0)
        //{
        //    particles->test[indx] = 1.0;
        //    particles->test0[indx] = particles->x[indx];
        //    particles->test1[indx] = particles->y[indx];
        //    particles->test2[indx] = particles->z[indx];
        //    particles->test3[indx] = signPartDotNormal;
        //}
        //else
        //{
        //    particles->test[indx] = particles->test[indx] + 1.0;
        //}
      } else {
        particles->hitWall[indx] = 2.0;
      }
    }
  }
};
#endif
