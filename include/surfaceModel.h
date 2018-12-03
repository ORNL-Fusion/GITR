#ifndef _SURFACE_
#define _SURFACE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particle.h"
#include <cmath>
#include <math.h>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#include <stdlib.h>
#endif

CUDA_CALLABLE_MEMBER
double screeningLength ( double Zprojectile, double Ztarget ) {
	double bohrRadius = 5.29177e-11;
	double screenLength;

	screenLength = 0.885341*bohrRadius*powf(powf(Zprojectile,(2/3)) + powf(Ztarget,(2/3)),(-1/2));

	return screenLength;
}

CUDA_CALLABLE_MEMBER
double stoppingPower (Particle p, double Mtarget, double Ztarget, double screenLength) {
	        double E0;
                double Q = 1.60217662e-19;
		double ke2 = 14.4e-10;
		double reducedEnergy;
	double stoppingPower;

	E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
	reducedEnergy = E0*(Mtarget/(p.amu+Mtarget))*(screenLength/(p.Z*Ztarget*ke2));
	stoppingPower = 0.5*log(1.0 + 1.2288*reducedEnergy)/(reducedEnergy + 0.1728*sqrt(reducedEnergy) + 0.008*pow(reducedEnergy, 0.1504));

	return stoppingPower;	
}

struct erosion { 

    const double dt;

    erosion(double _dt) : dt(_dt) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(Particle &p) const { 
	double screenLength;
	double stopPower;
	double q = 18.6006;
	double lambda = 2.2697;
    float eps =0.0; 
    double mu = 3.1273;
	double Eth = 24.9885;
	double Y0;
	double Ztarget = 74.0;
	double Mtarget = 183.84;
	double term;
	double E0;

	screenLength = screeningLength(p.Z, Ztarget);
	stopPower = stoppingPower(p, Mtarget, Ztarget, screenLength); 
	E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
	term = pow((E0/Eth - 1),mu);
	Y0 = q*stopPower*term/(lambda + term);
    	}
     
};

struct reflection {

    const double dt;
    int nLines;
    Boundary * boundaryVector;

    reflection(double _dt, int _nLines,Boundary * _boundaryVector) : dt(_dt),nLines(_nLines), boundaryVector(_boundaryVector) {}

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(Particle &p) const {
            if(p.hitWall == 1.0)
        	{
			float reducedEnergyMultiplier = 5e-7;//for W on W
            float E0 = 0.0;
            float reducedEnergy = 0.0;
			float a1 = -3.685;
			float a2 = 0.0278;
			float a3 = 7.825e-5;
			float a4 = -1.109;
			float Rn = 0.0;
			float Re = 0.0;
            float thetaImpact = 0.0;
            float thetaAzimuthal = 0.0;
            float particleTrackVector[3] = {0.0,0.0,0.0};
            float surfaceNormalVector[3] = {0.0,0.0,0.0};
            float surfaceNormalVectorRotated[3] = {0.0,0.0,0.0};
            float surfaceNormalVectorIn[3] = {0.0,0.0,0.0};
            float norm_part = 0.0;
            float norm_normal = 0.0;
            float partDotNormal = 0.0;
            float signPartDotNormal = 0.0;
		
		        E0 = 0.5*p.amu*1.6737236e-27*(p.vx*p.vx + p.vy*p.vy+ p.vz*p.vz)/1.60217662e-19;
			reducedEnergy = E0*reducedEnergyMultiplier;
            particleTrackVector[0] = p.vx;
            particleTrackVector[1] = p.vy;
            particleTrackVector[2] = p.vz;
            //surfaceNormalVector[0] = -boundaryVector[p.wallIndex].slope_dzdx;
            surfaceNormalVector[1] = 0.0;
            surfaceNormalVector[2] = 1.0;
            //std::cout << "velocities " << p.vx << " " << p.vy << " " << p.vz << std::endl;
            //std::cout << "surface norm " << surfaceNormalVector[0] << " " << surfaceNormalVector[1] << " " << surfaceNormalVector[2] << std::endl;
#if USECYLSYMM > 0
            thetaAzimuthal = atan2(p.y,p.x);   
            surfaceNormalVectorRotated[0] = cos(thetaAzimuthal)*surfaceNormalVector[0];
            surfaceNormalVectorRotated[1] = sin(thetaAzimuthal)*surfaceNormalVector[0];
#else
            surfaceNormalVectorRotated[0] = surfaceNormalVector[0];
            surfaceNormalVectorRotated[1] = surfaceNormalVector[1];
#endif
            surfaceNormalVectorRotated[2] = surfaceNormalVector[2];

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
    int surfaceHit = boundaryVector[wallHit].surfaceNumber;
    int surface = boundaryVector[wallHit].surface;
    float eInterpVal=0.0;
    float aInterpVal=0.0;
    float weight = particles->weight[indx];
    float vx = particles->vx[indx];
    float vy = particles->vy[indx];
    float vz = particles->vz[indx];
    #if FLUX_EA > 0
      float dEdist = (Edist - E0dist)/nEdist;
      float dAdist = (Adist - A0dist)/nAdist;
      int AdistInd=0;
      int EdistInd=0;
    #endif
        particles->firstCollision[indx]=1;
    particleTrackVector[0] = vx;
    particleTrackVector[1] = vy;
    particleTrackVector[2] = vz;
    norm_part = sqrt(particleTrackVector[0]*particleTrackVector[0] + particleTrackVector[1]*particleTrackVector[1] + particleTrackVector[2]*particleTrackVector[2]);
    E0 = 0.5*particles->amu[indx]*1.6737236e-27*(norm_part*norm_part)/1.60217662e-19;
    //std::cout << "Particle hit wall with energy " << E0 << std::endl;
    //std::cout << "Particle hit wall with v " << vx << " " << vy << " " << vz<< std::endl;
    //std::cout << "Particle amu norm_part " << particles->amu[indx] << " " << vy << " " << vz<< std::endl;
    wallIndex = particles->wallIndex[indx];
    boundaryVector[wallHit].getSurfaceNormal(surfaceNormalVector);
            particleTrackVector[0] = particleTrackVector[0]/norm_part;
            particleTrackVector[1] = particleTrackVector[1]/norm_part;
            particleTrackVector[2] = particleTrackVector[2]/norm_part;

            partDotNormal = particleTrackVector[0]*surfaceNormalVectorRotated[0] + particleTrackVector[1]*surfaceNormalVectorRotated[1] + particleTrackVector[2]*surfaceNormalVectorRotated[2];
            thetaImpact = acos(partDotNormal);
            signPartDotNormal = sgn(partDotNormal);
            surfaceNormalVectorRotated[0] = -surfaceNormalVectorRotated[0]*signPartDotNormal;
            surfaceNormalVectorRotated[1] = -surfaceNormalVectorRotated[1]*signPartDotNormal;
            surfaceNormalVectorRotated[2] = -surfaceNormalVectorRotated[2]*signPartDotNormal;

            Rn = exp(a1*pow(reducedEnergy,a2))/(1 + exp(a3*pow(reducedEnergy,a4)));
            //std::cout << "calculated Rn " << Rn << std::endl;
            Rn = 0.25;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
                curandState tmpState;
                float r7 = curand_uniform(&p.streams[6]);
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
                  aInterpVal = interp3d (r8,thetaImpact,log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                    angleDistGrid01,A_sputtRefDistIn,
                                    E_sputtRefDistIn,ADist_CDF_R_regrid);
                   eInterpVal = interp3d ( r9,thetaImpact,log10(E0),
                           nE_sputtRefDistOutRef,nA_sputtRefDistIn,nE_sputtRefDistIn,
                                         energyDistGrid01Ref,A_sputtRefDistIn,
                                         E_sputtRefDistIn,EDist_CDF_R_regrid );
                   //newWeight=(R0/(1.0f-sputtProb))*weight;
		   newWeight = weight*(totalYR);
    #if FLUX_EA > 0
              EdistInd = floor((eInterpVal-E0dist)/dEdist);
              AdistInd = floor((aInterpVal-A0dist)/dAdist);
              if((EdistInd >= 0) && (EdistInd < nEdist) && 
                 (AdistInd >= 0) && (AdistInd < nAdist))
              {
                  surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                    surfaces->reflDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
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
                  aInterpVal = interp3d(r8,thetaImpact,log10(E0),
                          nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
                          angleDistGrid01,A_sputtRefDistIn,
                          E_sputtRefDistIn,ADist_CDF_Y_regrid);
                  eInterpVal = interp3d(r9,thetaImpact,log10(E0),
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
              EdistInd = floor((eInterpVal-E0dist)/dEdist);
              AdistInd = floor((aInterpVal-A0dist)/dAdist);
              if((EdistInd >= 0) && (EdistInd < nEdist) && 
                 (AdistInd >= 0) && (AdistInd < nAdist))
              {
                //std::cout << " particle sputters with " << EdistInd << AdistInd <<  std::endl;
                  surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] = 
                    surfaces->sputtDistribution[surfaceHit*nEdist*nAdist + EdistInd*nAdist + AdistInd] +  newWeight;
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
                EdistInd = floor((E0-E0dist)/dEdist);
                AdistInd = floor((thetaImpact-A0dist)/dAdist);
              
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
	    if( boundaryVector[wallHit].Z > 0.0 && newWeight > 0.0)
            {
                particles->weight[indx] = newWeight;
                particles->hitWall[indx] = 0.0;
                particles->charge[indx] = 0.0;
                float V0 = sqrt(2*eInterpVal*1.602e-19/(particles->amu[indx]*1.66e-27));
                particles->newVelocity[indx] = V0;
    vSampled[0] = V0*sin(aInterpVal*3.1415/180)*cos(2.0*3.1415*r10);
    vSampled[1] = V0*sin(aInterpVal*3.1415/180)*sin(2.0*3.1415*r10);
    vSampled[2] = V0*cos(aInterpVal*3.1415/180);
    boundaryVector[wallHit].transformToSurface(vSampled);
    particles->vx[indx] = -signPartDotNormal*vSampled[0];
    particles->vy[indx] = -signPartDotNormal*vSampled[1];
    particles->vz[indx] = -signPartDotNormal*vSampled[2];
            //if(particles->test[indx] == 0.0)
            //{
            //    particles->test[indx] = 1.0;
            //    particles->test0[indx] = aInterpVal;
            //    particles->test1[indx] = eInterpVal;
            //    particles->test2[indx] = V0;
            //    particles->test3[indx] = vSampled[2];
            //}            

    particles->xprevious[indx] = particles->x[indx] + -signPartDotNormal*surfaceNormalVector[0]*1e-5;
    particles->yprevious[indx] = particles->y[indx] + -signPartDotNormal*surfaceNormalVector[1]*1e-5;
    particles->zprevious[indx] = particles->z[indx] + -signPartDotNormal*surfaceNormalVector[2]*1e-5;
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
            }
                else
                {
                    //std::cout << " no reflection" << std::endl;
                    p.hitWall = 2.0;
                }
		}
	}
};
#endif
