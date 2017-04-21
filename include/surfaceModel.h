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
    Particles * particles;
    const double dt;
    int nLines;
    Boundary * boundaryVector;
    int nAngle;
    int nEnergy;
    float* spYlGridAngle;
    float* spYlGridE;
    float* spYl;

    reflection(Particles* _particles, double _dt, int _nLines,Boundary * _boundaryVector,
            int _nAngle, int _nEnergy, float* _spYlGridAngle, float* _spYlGridE, float* _spYl) : 
        particles(_particles), dt(_dt),nLines(_nLines), boundaryVector(_boundaryVector),
        nAngle(_nAngle),nEnergy(_nEnergy),spYlGridAngle(_spYlGridAngle),
        spYlGridE(_spYlGridE), spYl(_spYl) {}

CUDA_CALLABLE_MEMBER_DEVICE
void operator()(std::size_t indx) const {
            if(particles->hitWall[indx] == 1.0)
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
		
		        E0 = 0.5*particles->amu[indx]*1.6737236e-27*(particles->vx[indx]*particles->vx[indx] + particles->vy[indx]*particles->vy[indx]+ particles->vz[indx]*particles->vz[indx])/1.60217662e-19;
			reducedEnergy = E0*reducedEnergyMultiplier;
            particleTrackVector[0] = particles->vx[indx];
            particleTrackVector[1] = particles->vy[indx];
            particleTrackVector[2] = particles->vz[indx];
            //surfaceNormalVector[0] = -boundaryVector[particles->wallIndex].slope_dzdx;
            surfaceNormalVector[1] = 0.0;
            surfaceNormalVector[2] = 1.0;
            //std::cout << "velocities " << particles->vx[indx] << " " << particles->vy[indx] << " " << particles->vz[indx] << std::endl;
            //std::cout << "surface norm " << surfaceNormalVector[0] << " " << surfaceNormalVector[1] << " " << surfaceNormalVector[2] << std::endl;
#if USECYLSYMM > 0
            thetaAzimuthal = atan2(particles->y[indx],particles->x[indx]);   
            surfaceNormalVectorRotated[0] = cos(thetaAzimuthal)*surfaceNormalVector[0];
            surfaceNormalVectorRotated[1] = sin(thetaAzimuthal)*surfaceNormalVector[0];
#else
            surfaceNormalVectorRotated[0] = surfaceNormalVector[0];
            surfaceNormalVectorRotated[1] = surfaceNormalVector[1];
#endif
            surfaceNormalVectorRotated[2] = surfaceNormalVector[2];

            norm_part = sqrt(particleTrackVector[0]*particleTrackVector[0] + particleTrackVector[1]*particleTrackVector[1] + particleTrackVector[2]*particleTrackVector[2]);
            norm_normal = sqrt(surfaceNormalVectorRotated[0]*surfaceNormalVectorRotated[0] + surfaceNormalVectorRotated[1]*surfaceNormalVectorRotated[1] + surfaceNormalVectorRotated[2]*surfaceNormalVectorRotated[2]);
    //        std::cout << "norm of particle track " << norm_part << std::endl;
    //        std::cout << "norm of surface normal " << norm_normal << std::endl;

            surfaceNormalVectorRotated[0] = surfaceNormalVectorRotated[0]/norm_normal;
            surfaceNormalVectorRotated[1] = surfaceNormalVectorRotated[1]/norm_normal;
            surfaceNormalVectorRotated[2] = surfaceNormalVectorRotated[2]/norm_normal;
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
                float r7 = curand_uniform(&particles->streams[6]);
#else
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                float r7 = dist(particles->streams[6]);
#endif

#else
                float r7 = 0.0;
#endif
		        if(r7 <= Rn)
        		{
      //              std::cout << "particle reflected " << std::endl;
		//		    std::cout << "from surface " << particles->wallIndex << " material " <<boundaryVector[particles->wallIndex].Z << std::endl;
                float a1_e = -7.168;
                float a2_e = 0.01685;
                float a3_e = 7.005e-5;
                float a4_e = -1.343;
				float Re = exp(a1_e*pow(reducedEnergy,a2_e))/(1 + exp(a3_e*pow(reducedEnergy,a4_e))); 
				float launch_unitVector[3] = {0, 0, 1.0};
				float launchEnergy = E0*Re/Rn;
                launchEnergy = 10.0;
              //  std::cout << "Re, E0 and launch energy" << Re << " " << E0 << " " << launchEnergy << std::endl;
                particles->charge[indx] = 0.0;
				particles->vx[indx] = sqrt(2*launchEnergy*1.60217662e-19/(particles->amu[indx]*1.6737236e-27))*surfaceNormalVectorRotated[0];//launchEnergy*launch_unitVector[0];
				particles->vy[indx] = sqrt(2*launchEnergy*1.60217662e-19/(particles->amu[indx]*1.6737236e-27))*surfaceNormalVectorRotated[1];//launchEnergy*launch_unitVector[1];
				particles->vz[indx] = sqrt(2*launchEnergy*1.60217662e-19/(particles->amu[indx]*1.6737236e-27))*surfaceNormalVectorRotated[2];
                //particles->x = 1.33;//particles->xprevious[indx];
                //particles->y = 0.0;
                //particles->z = 1.34;//boundaryVector[particles->wallIndex].z1 + 0.01;
                particles->xprevious[indx] = particles->xprevious[indx] + 0.0001*surfaceNormalVectorRotated[0];
                particles->yprevious[indx] = particles->yprevious[indx] + 0.0001*surfaceNormalVectorRotated[1];
                particles->zprevious[indx] = particles->zprevious[indx] + 0.0001*surfaceNormalVectorRotated[2];//boundaryVector[particles->wallIndex].z1 + 0.01;
                //std::cout << "particle launched from " << particles->xprevious[indx] << particles->yprevious[indx] << particles->zprevious[indx] << std::endl;
                //std::cout << "particle launched at velocity " << particles->vx[indx] << particles->vy[indx] << particles->vz[indx] << std::endl;

				particles->hitWall[indx] = 0.0;
        		}
                else
                {
                    //std::cout << " no reflection" << std::endl;
                    particles->hitWall[indx] = 2.0;
                }
		}
	}
};
#endif
