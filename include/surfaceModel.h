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

	E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/1.60217662e-19;
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
    float eps =0.0; 
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
    Particles * particlesArray;
    const double dt;
    int nLines;
    Boundary * boundaryVector;
    int nAngle;
    int nEnergy;
    float* spYlGridAngle;
    float* spYlGridE;
    float* spYl;
    int nSegmentsAngle;
    float* sourceAngleSegments;
    float* angleCDF;
    int nThompDistPoints;
    float max_Energy;
    float* CDFThompson;

    reflection(Particles* _particlesArray, double _dt, int _nLines,Boundary * _boundaryVector,
            int _nAngle, int _nEnergy, float* _spYlGridAngle, float* _spYlGridE, float* _spYl,
            int _nSegmentsAngle, float* _sourceAngleSegments, float* _angleCDF,
            int _nThompDistPoints, float _max_Energy, float* _CDFThompson) : 
        particlesArray(_particlesArray), dt(_dt),nLines(_nLines), boundaryVector(_boundaryVector),
        nAngle(_nAngle),nEnergy(_nEnergy),spYlGridAngle(_spYlGridAngle),
        spYlGridE(_spYlGridE), spYl(_spYl),
      nSegmentsAngle(_nSegmentsAngle), sourceAngleSegments(_sourceAngleSegments),
   angleCDF(_angleCDF), 
   nThompDistPoints(_nThompDistPoints),max_Energy(_max_Energy), CDFThompson(_CDFThompson) {}

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const {
    
            if(particlesArray->hitWall[indx] == 1.0)
            {
                particlesArray->hitWall[indx] = 0.0;
            }
            //if(particles->hitWall[indx] == 1.0)
             //      if( boundaryVector[particles->wallHit[indx]].Z > 0.0)
              //     {
        	//{
                
	        /*
                //float reducedEnergyMultiplier = 5e-7;//for W on W
            float E0 = 0.0;
            float thetaImpact = 0.0;
            //float thetaAzimuthal = 0.0;
            float particleTrackVector[3] = {0.0,0.0,0.0};
            float surfaceNormalVector[3] = {0.0,0.0,0.0};
            float norm_part = 0.0;
            float norm_normal = 0.0;
            
            float partDotNormal = 0.0;
            float signPartDotNormal = 0.0;
            float Enew = 0.0f;
            float angleSample = 0.0f;
            int wallIndex = 0;

		        E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/1.60217662e-19;
		//	reducedEnergy = E0*reducedEnergyMultiplier;
            particleTrackVector[0] = particles->vx[indx];
            particleTrackVector[1] = particles->vy[indx];
            particleTrackVector[2] = particles->vz[indx];
            //surfaceNormalVector[0] = -boundaryVector[particles->wallIndex].slope_dzdx;
            wallIndex = particles->wallIndex[indx];
            norm_normal = boundaryVector[wallIndex].plane_norm; 
            surfaceNormalVector[0] = boundaryVector[wallIndex].a/norm_normal;
            surfaceNormalVector[1] = boundaryVector[wallIndex].b/norm_normal;
            
            surfaceNormalVector[2] = boundaryVector[wallIndex].c/norm_normal;
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
            signPartDotNormal = sgn(partDotNormal);
            float Y0 = interp2dUnstructured(thetaImpact*180.0/3.1415,E0,nAngle,nEnergy, 
                    spYlGridAngle,spYlGridE, spYl);
            
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
                float r7 = curand_uniform(&particles->streams_surf[indx]);
                float r8 = curand_uniform(&particles->streams_surf[indx]);
#else
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                float r7 = dist(particles->streams_surf[indx]);
                float r8 = dist(particles->streams_surf[indx]);
#endif

#else
#if __CUDACC__
                    float r1 = curand_uniform(state);
#else
                            std::uniform_real_distribution<float> dist(0.0, 1.0);
                                        float r1=dist(state[0]);
#endif
                //float r7 = 0.0;
#endif

                //deposit on surface
                int wallHit = particles->wallHit[indx];
#if USESURFACEMODEL == 0 
#if USE_CUDA > 0
                atomicAdd(&boundaryVector[wallHit].impacts, particles->weight[indx]);
#else
               boundaryVector[wallHit].impacts = boundaryVector[wallHit].impacts +  particlesPointer->weight[indx];
#endif
#endif 

                //reflect with weight and new initial conditions
                particles->weight[indx] = particles->weight[indx]*Y0;
*/
               // particles->hitWall[indx] = 0.0;
/*
                Enew = interp1dUnstructured(r7,nThompDistPoints, max_Energy, CDFThompson);
                
                angleSample = interp1dUnstructured2(r8,nSegmentsAngle,sourceAngleSegments , angleCDF);
                 float V0 = sqrt(2*E0*1.602e-19/(particles->amu[indx]*1.66e-27));
    particles->vx[indx] = -signPartDotNormal*V0*surfaceNormalVector[0];
    particles->vy[indx] = -signPartDotNormal*V0*surfaceNormalVector[1];
    particles->vz[indx] = -signPartDotNormal*V0*surfaceNormalVector[2];
*/
            //}
              //     }
}
};
#endif
