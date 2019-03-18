#ifndef _PARDIFFUSION_
  #define _PARDIFFUSION_

  #ifdef __CUDACC__
    #define CUDA_CALLABLE_MEMBER __host__ __device__
    #define CUDA_CALLABLE_MEMBER_DEVICE __device__
  #else
    #define CUDA_CALLABLE_MEMBER
    #define CUDA_CALLABLE_MEMBER_DEVICE
  #endif
  
  #include "Particles.h"
  #include "math.h"
  #ifdef __CUDACC__
    #include <thrust/random.h>
  #else
    #include <random>
    #include <stdlib.h>
  #endif

  struct parDiffusion { 
    Particles *particlesPointer;
    const float dt;
    int nR_flowV;
    int nY_flowV;
    int nZ_flowV;
    float* flowVGridr;
    float* flowVGridy;
    float* flowVGridz;
    float* flowVr;
    float* flowVz;
    float* flowVt;
    int nR_Dens;
    int nZ_Dens;
    float* DensGridr;
    float* DensGridz;
    float* ni;
    int nR_Temp;
    int nZ_Temp;
    float* TempGridr;
    float* TempGridz;
    float* ti;
    float* te;
    float background_Z;
    float background_amu;
    int nR_Bfield;
    int nZ_Bfield;
    float * BfieldGridR;
    float * BfieldGridZ;
    float * BfieldR;
    float * BfieldZ;
    float * BfieldT;
    float dv[3];
    #if __CUDACC__
      curandState *state;
    #else
      std::mt19937 *state;
    #endif

    parDiffusion(Particles *_particlesPointer,float _dt, 
    #if __CUDACC__
      curandState *_state,
    #else
      std::mt19937 *_state,
    #endif
    int _nR_flowV,int _nY_flowV, int _nZ_flowV,    
    float* _flowVGridr,float* _flowVGridy,
    float* _flowVGridz,float* _flowVr,
    float* _flowVz,float* _flowVt,
    int _nR_Dens,int _nZ_Dens,float* _DensGridr,
    float* _DensGridz,float* _ni,int _nR_Temp, int _nZ_Temp,
    float* _TempGridr, float* _TempGridz,float* _ti,float* _te,
    float _background_Z, float _background_amu,
    int _nR_Bfield, int _nZ_Bfield,
    float * _BfieldGridR ,float * _BfieldGridZ ,
    float * _BfieldR ,float * _BfieldZ ,
    float * _BfieldT )
    : particlesPointer(_particlesPointer), dt(_dt),state(_state), 
    nR_flowV(_nR_flowV),nY_flowV(_nY_flowV), nZ_flowV(_nZ_flowV), 
    flowVGridr(_flowVGridr),flowVGridy(_flowVGridy),
    flowVGridz(_flowVGridz), flowVr(_flowVr),flowVz(_flowVz), flowVt(_flowVt),
    nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), 
    DensGridz(_DensGridz),ni(_ni),
    nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), 
    TempGridz(_TempGridz),
    ti(_ti),te(_te),background_Z(_background_Z), 
    background_amu(_background_amu),
    nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridR(_BfieldGridR), 
    BfieldGridZ(_BfieldGridZ),BfieldR(_BfieldR), BfieldZ(_BfieldZ), 
    BfieldT(_BfieldT),
    dv{0.0f,0.0f,0.0f} {} 

    CUDA_CALLABLE_MEMBER_DEVICE    
    void operator()(std::size_t indx)  { 
      if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0)
      { 
        float pi = 3.14159265;   
	float k_boltz = 1.38e-23*11604/1.66e-27;
	float T_background = 0.0;
        float nu_friction = 0.0;
        float nu_deflection = 0.0;
        float nu_parallel = 0.0;
        float nu_energy = 0.0;
        float flowVelocity[3]= {0.0f};
        float vUpdate[3]= {0.0f};
        float relativeVelocity[3] = {0.0f};
        float velocityCollisions[3]= {0.0f};	
        float velocityRelativeNorm;	
        float parallel_direction[3] = {0.0f};
        float perp_direction1[3] = {0.0f};
        float perp_direction2[3] = {0.0f};
        float parallel_contribution;
        float dv_perp1[3] = {0.0f};
        float dv_perp2[3] = {0.0f};
        float x = particlesPointer->xprevious[indx];
        float y = particlesPointer->yprevious[indx];
        float z = particlesPointer->zprevious[indx];
        float vx = particlesPointer->vx[indx];
        float vy = particlesPointer->vy[indx];
        float vz = particlesPointer->vz[indx];
        float vPartNorm = 0.0f;
	float Bnorm = 0.0f;
        #if FLOWV_INTERP == 3 
          interp3dVector (&flowVelocity[0], particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_flowV,nY_flowV,nZ_flowV,
                flowVGridr,flowVGridy,flowVGridz,flowVr,flowVz,flowVt);
        #elif FLOWV_INTERP < 3    
          #if USEFIELDALIGNEDVALUES > 0
            interpFieldAlignedVector(&flowVelocity[0],
                                 particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                                 nR_flowV,nZ_flowV,
                                 flowVGridr,flowVGridz,flowVr,
                                 flowVz,flowVt,nR_Bfield,nZ_Bfield,
                                 BfieldGridR,BfieldGridZ,BfieldR,
                                 BfieldZ,BfieldT);
          #else
            interp2dVector(flowVelocity,particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                        nR_flowV,nZ_flowV,
                        flowVGridr,flowVGridz,flowVr,flowVz,flowVt);
          #endif
        #endif
        vPartNorm = sqrt(vx*vx + vy*vy + vz*vz);
        relativeVelocity[0] = vx - flowVelocity[0];
        relativeVelocity[1] = vy - flowVelocity[1];
        relativeVelocity[2] = vz - flowVelocity[2];
        float vRel2 = relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2];
        velocityRelativeNorm = vectorNorm(relativeVelocity);

        #if PARTICLESEEDS > 0
          #ifdef __CUDACC__
            float r1 = (curand_uniform(&state[indx]) - 0.5);
          #else
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            float r1 = 2.0*floor(dist(state[indx]) + 0.5)-1.0;
          #endif
          #else
          #if __CUDACC__
          #else
            std::uniform_real_distribution<float> dist(0.0, 1.0);
          #endif
        #endif 
	getSlowDownFrequencies (nu_friction,nu_deflection,nu_parallel, nu_energy,
                                x,y,z,
                                vx,vy,vz,
                                particlesPointer->charge[indx],particlesPointer->amu[indx],
                                nR_flowV,  nZ_flowV, flowVGridr,
                                flowVGridz,flowVr,
                                flowVz,flowVt,
                                nR_Dens, nZ_Dens,DensGridr,
                                DensGridz,ni, nR_Temp,  nZ_Temp,
                                TempGridr, TempGridz,ti,te, background_Z, background_amu,
                                nR_Bfield,
                                nZ_Bfield,
                                BfieldGridR,
                                BfieldGridZ,
                                BfieldR,
                                BfieldZ,
                                BfieldT,T_background);

            interp2dVector(parallel_direction,particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                        nR_Bfield,nZ_Bfield,
                        BfieldGridR,BfieldGridZ,BfieldR,BfieldZ,BfieldT);
	    float Dpar = vPartNorm*vPartNorm/nu_friction;
	    //std::cout << "vPartNorm " << vPartNorm << endl;
	    //std::cout << "nu_friction " << nu_friction << endl;
	    //std::cout << "Dpar " << Dpar << endl;
	    vectorNormalize(parallel_direction,parallel_direction);
            if(vectorNorm(parallel_direction) == 0.0){
	      Dpar=0.0;
	      parallel_direction[0]=0.0;
	      parallel_direction[1]=0.0;
	      parallel_direction[2]=0.0;
	    }
            if(nu_friction == 0.0){Dpar=0.0;}
	    particlesPointer->xprevious[indx] = particlesPointer->xprevious[indx] + r1*sqrt(Dpar*dt)*parallel_direction[0];
	    particlesPointer->yprevious[indx] = particlesPointer->yprevious[indx] + r1*sqrt(Dpar*dt)*parallel_direction[1];
	    particlesPointer->zprevious[indx] = particlesPointer->zprevious[indx] + r1*sqrt(Dpar*dt)*parallel_direction[2];

	}
    }
     
};

#endif
