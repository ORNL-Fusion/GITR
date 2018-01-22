#ifndef _CFDIFFUSION_
#define _CFDIFFUSION_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include <cmath>

struct crossFieldDiffusion { 
    Particles *particlesPointer;
    const float dt;
	const float diffusionCoefficient;
    int nR_Bfield;
    int nZ_Bfield;
    float * BfieldGridRDevicePointer;
    float * BfieldGridZDevicePointer;
    float * BfieldRDevicePointer;
    float * BfieldZDevicePointer;
    float * BfieldTDevicePointer;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    crossFieldDiffusion(Particles *_particlesPointer, float _dt,
#if __CUDACC__
                            curandState *_state,
#else
                                            std::mt19937 *_state,
#endif
            float _diffusionCoefficient,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridRDevicePointer,float * _BfieldGridZDevicePointer,
            float * _BfieldRDevicePointer,float * _BfieldZDevicePointer,
            float * _BfieldTDevicePointer)
        : particlesPointer(_particlesPointer), dt(_dt),state(_state), diffusionCoefficient(_diffusionCoefficient),nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridRDevicePointer(_BfieldGridRDevicePointer), BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
       BfieldRDevicePointer(_BfieldRDevicePointer), BfieldZDevicePointer(_BfieldZDevicePointer), BfieldTDevicePointer(_BfieldTDevicePointer) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 

	    if(particlesPointer->hitWall[indx] == 0.0)
        {
           if(particlesPointer->charge[indx] > 0.0)
           { 
       
	        float perpVector[3]= {0, 0, 0};
	        float B[3] = {0.0,0.0,0.0};
            float Bmag = 0.0;
		float B_unit[3] = {0.0, 0.0, 0.0};
		float phi_random;
		float norm;
		float step;
        interp2dVector(&B[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_Bfield,nZ_Bfield,
                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
        Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
        	float r3 = curand_uniform(&state[indx]);
#else
        	std::uniform_real_distribution<float> dist(0.0, 1.0);
        	float r3=dist(state[indx]);
#endif 
#else
#if __CUDACC__
            float r3 = curand_uniform(&state[2]);
#else
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            float r3=dist(state[2]);
#endif
#endif
		phi_random = 2*3.14159265*r3;
		perpVector[0] = cos(phi_random);
		perpVector[1] = sin(phi_random);
		perpVector[2] = (-perpVector[0]*B_unit[0] - perpVector[1]*B_unit[1])/B_unit[2];

		if (B_unit[2] == 0){
			perpVector[2] = perpVector[1];
			perpVector[1] = (-perpVector[0]*B_unit[0] - perpVector[2]*B_unit[2])/B_unit[1];
		}
		
		if ((B_unit[0] == 1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0) || (B_unit[0] == -1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0))
		{
			perpVector[2] = perpVector[0];
			perpVector[0] = 0;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==1.0 && B_unit[2] ==0.0) || (B_unit[0] == 0.0 && B_unit[1] ==-1.0 && B_unit[2] ==0.0))
		{
			perpVector[1] = 0.0;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==1.0) || (B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==-1.0))
		{
			perpVector[2] = 0;
		}
		
		norm = sqrt(perpVector[0]*perpVector[0] + perpVector[1]*perpVector[1] + perpVector[2]*perpVector[2]);
		perpVector[0] = perpVector[0]/norm;
		perpVector[1] = perpVector[1]/norm;
		perpVector[2] = perpVector[2]/norm;
		
		step = sqrt(6*diffusionCoefficient*dt);

		particlesPointer->x[indx] = particlesPointer->xprevious[indx] + step*perpVector[0];
		particlesPointer->y[indx] = particlesPointer->yprevious[indx] + step*perpVector[1];
		particlesPointer->z[indx] = particlesPointer->zprevious[indx] + step*perpVector[2];
    	}
    } }
};

#endif
