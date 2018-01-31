#ifndef _COULOMB_
#define _COULOMB_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include <cmath>
#include <math.h>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#include <stdlib.h>
#endif

CUDA_CALLABLE_MEMBER
void getSlowDownFrequencies ( float& nu_friction, float& nu_deflection, float& nu_parallel,
			 	float& nu_energy, float x, float y,float z, float vx, float vy, float vz,float charge, float amu, 
                
    int nR_flowV,
    int nZ_flowV,
    float* flowVGridr,
    float* flowVGridz,
    float* flowVr,
    float* flowVz,
    float* flowVt,
    int nR_Dens,
    int nZ_Dens,
    float* DensGridr,
    float* DensGridz,
    float* ni,
    int nR_Temp,
    int nZ_Temp,
    float* TempGridr,
    float* TempGridz,
    float* ti, float* te,float background_Z, float background_amu
                ) {
        float Q = 1.60217662e-19;
        float EPS0 = 8.854187e-12;
	float pi = 3.14159265;
float MI = 1.6737236e-27;	
float ME = 9.10938356e-31;
        float te_eV = interp2dCombined(x,y,z,nR_Temp,nZ_Temp,TempGridr,TempGridz,te);
        float ti_eV = interp2dCombined(x,y,z,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti);
            float density = interp2dCombined(x,y,z,nR_Dens,nZ_Dens,DensGridr,DensGridz,ni);
            //std::cout << "ion t and n " << te_eV << "  " << density << std::endl;
    float flowVelocity[3]= {0, 0, 0};
	float relativeVelocity[3] = {0.0, 0.0, 0.0};
	float velocityNorm = 0.0f;
	float lam_d;
	float lam;
	float gam_electron_background;
	float gam_ion_background;
	float a_electron = 0.0;
	float a_ion = 0.0;
	float xx;
	float psi_prime;
	float psi_psiprime;
	float psi;
	float xx_e;
	float psi_prime_e;
	float psi_psiprime_e;
	float psi_e;
	float nu_0_i;
	float nu_0_e;
    float nu_friction_i; 
    float nu_deflection_i;
    float nu_parallel_i; 
    float nu_energy_i;
    float nu_friction_e; 
    float nu_deflection_e;
    float nu_parallel_e; 
    float nu_energy_e;
                interp2dVector(&flowVelocity[0],x,y,z,nR_flowV,nZ_flowV,
                                       flowVGridr,flowVGridz,flowVr,flowVz,flowVt);
	relativeVelocity[0] = vx - flowVelocity[0];
	relativeVelocity[1] = vy - flowVelocity[1];
	relativeVelocity[2] = vz - flowVelocity[2];
	velocityNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);                
	    //std::cout << "velocity norm " << velocityNorm << std::endl;	
    //for(int i=1; i < nSpecies; i++)
		//{
			lam_d = sqrtf(EPS0*te_eV/(density*powf(background_Z,2)*Q));//only one q in order to convert to J
                	lam = 4.0*pi*density*powf(lam_d,3);
                	gam_electron_background = 0.238762895*powf(charge,2)*logf(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)
                	gam_ion_background = 0.238762895*powf(charge,2)*powf(background_Z,2)*logf(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)
                    //std::cout << "gam components " <<gam_electron_background << " " << pow(Q,4) << " " << " " << pow(background_Z,2) << " " << log(lam)<< std::endl; 
                	a_ion = background_amu*MI/(2*ti_eV*Q);// %q is just to convert units - no z needed
                	a_electron = ME/(2*te_eV*Q);// %q is just to convert units - no z needed
                
                	xx = powf(velocityNorm,2)*a_ion;
                	psi_prime = 2*sqrtf(xx/pi)*expf(-xx);
                	psi_psiprime = erf(1.0*sqrtf(xx));
                	psi = psi_psiprime - psi_prime;
                	xx_e = powf(velocityNorm,2)*a_electron;
                	psi_prime_e = 2*sqrtf(xx_e/pi)*expf(-xx_e);
                	psi_psiprime_e = erf(1.0*sqrtf(xx_e));
                	psi_e = psi_psiprime_e - psi_prime_e;
                	nu_0_i = gam_electron_background*density/powf(velocityNorm,3);
                	nu_0_e = gam_ion_background*density/powf(velocityNorm,3);
                	nu_friction_i = (1+amu/background_amu)*psi*nu_0_i;
                	nu_deflection_i = 2*(psi_psiprime - psi/(2*xx))*nu_0_i;
                	nu_parallel_i = psi/xx*nu_0_i;
                	nu_energy_i = 2*(amu/background_amu*psi - psi_prime)*nu_0_i;
                	nu_friction_e = (1+amu/(ME/MI))*psi_e*nu_0_e;
                	nu_deflection_e = 2*(psi_psiprime_e - psi_e/(2*xx_e))*nu_0_e;
                	nu_parallel_e = psi_e/xx_e*nu_0_e;
                	nu_energy_e = 2*(amu/(ME/MI)*psi_e - psi_prime_e)*nu_0_e;

                    //std::cout << "lam_d lam gami game ai ae" << lam_d << " " << lam << " " << gam_ion_background << " " << gam_electron_background << " " << a_ion << " " << a_electron << std::endl;
                    //std::cout << "x psi_prime psi_psiprime psi" << x << " " << psi_prime << " " << psi_psiprime << " " << psi << " " << nu_0_i << std::endl;
                    //std::cout << "nu friction, parallel perp energy " << nu_friction << " " << nu_parallel << " " <<nu_deflection << " " << nu_energy << std::endl;
	//	}
    nu_friction = nu_friction_i + nu_friction_e;
    nu_deflection = nu_deflection_i + nu_deflection_e;
    nu_parallel = nu_parallel_i + nu_parallel_e;
    nu_energy = nu_energy_i + nu_energy_e;

    //std::cout << "nu friction i e total " << nu_deflection_i << " " << nu_deflection_e << " " <<nu_deflection  << std::endl;

}

CUDA_CALLABLE_MEMBER
void getSlowDownDirections (float parallel_direction[], float perp_direction1[], float perp_direction2[],
        float xprevious, float yprevious, float zprevious,float vx, float vy, float vz,
    int nR_flowV,
    int nY_flowV,
    int nZ_flowV,
    float* flowVGridr,
    float* flowVGridy,
    float* flowVGridz,
    float* flowVr,
    float* flowVz,
    float* flowVt,
    
                        int nR_Bfield, int nZ_Bfield,
                        float* BfieldGridR ,float* BfieldGridZ ,
                        float* BfieldR ,float* BfieldZ ,
                 float* BfieldT 
    ) {
	        float flowVelocity[3]= {0, 0, 0.0};
                float relativeVelocity[3] = {0.0, 0.0, 0.0};
                float B[3] = {};
                float Bmag = 0.0;
		float B_unit[3] = {0.0, 0.0, -1.0};
		float velocityRelativeNorm;
		float s1;
		float s2;
                interp2dVector(&B[0],xprevious,yprevious,zprevious,nR_Bfield,nZ_Bfield,
                                       BfieldGridR,BfieldGridZ,BfieldR,BfieldZ,BfieldT);
        Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;

#if USE3DTETGEOM
        interp3dVector (&flowVelocity[0], xprevious,yprevious,zprevious,nR_flowV,nY_flowV,nZ_flowV,
                flowVGridr,flowVGridy,flowVGridz,flowVr,flowVz,flowVt);
#else    
                interp2dVector(&flowVelocity[0],xprevious,yprevious,zprevious,nR_flowV,nZ_flowV,
                                       flowVGridr,flowVGridz,flowVr,flowVz,flowVt);
#endif
                relativeVelocity[0] = vx - flowVelocity[0];
                relativeVelocity[1] = vy - flowVelocity[1];
                relativeVelocity[2] = vz - flowVelocity[2];
                velocityRelativeNorm = sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);

		parallel_direction[0] = relativeVelocity[0]/velocityRelativeNorm;
		parallel_direction[1] = relativeVelocity[1]/velocityRelativeNorm;
		parallel_direction[2] = relativeVelocity[2]/velocityRelativeNorm;

		s1 = parallel_direction[0]*B_unit[0]+parallel_direction[1]*B_unit[1]+parallel_direction[2]*B_unit[2];
            	s2 = sqrt(1-s1*s1);
            
            	perp_direction1[0] = 1/s2*(s1*parallel_direction[0] - B_unit[0]);
		perp_direction1[1] = 1/s2*(s1*parallel_direction[1] - B_unit[1]);
		perp_direction1[2] = 1/s2*(s1*parallel_direction[2] - B_unit[2]);
           
                perp_direction2[0] = 1/s2*(parallel_direction[1]*B_unit[2] - parallel_direction[2]*B_unit[1]);
                perp_direction2[1] = 1/s2*(parallel_direction[2]*B_unit[0] - parallel_direction[0]*B_unit[2]);
                perp_direction2[2] = 1/s2*(parallel_direction[0]*B_unit[1] - parallel_direction[1]*B_unit[0]); 
            if (s2 == 0)
		{
                perp_direction2[0] = parallel_direction[2];
		perp_direction2[1] = parallel_direction[0];
		perp_direction2[2] =  parallel_direction[1];

                s1 = parallel_direction[0]*perp_direction2[0]+parallel_direction[1]*perp_direction2[1]+parallel_direction[2]*perp_direction2[2];
                s2 = sqrt(1-s1*s1);
		perp_direction1[0] = -1/s2*(parallel_direction[1]*perp_direction2[2] - parallel_direction[2]*perp_direction2[1]);
                perp_direction1[1] = -1/s2*(parallel_direction[2]*perp_direction2[0] - parallel_direction[0]*perp_direction2[2]);
                perp_direction1[2] = -1/s2*(parallel_direction[0]*perp_direction2[1] - parallel_direction[1]*perp_direction2[0]);
            }
	
}

struct coulombCollisions { 
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
#if __CUDACC__
            curandState *state;
#else
            std::mt19937 *state;
#endif

    coulombCollisions(Particles *_particlesPointer,float _dt, 
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nR_flowV,int _nY_flowV, int _nZ_flowV,    float* _flowVGridr,float* _flowVGridy,
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
        : particlesPointer(_particlesPointer), dt(_dt),state(_state), nR_flowV(_nR_flowV),nY_flowV(_nY_flowV), nZ_flowV(_nZ_flowV), flowVGridr(_flowVGridr),flowVGridy(_flowVGridy),
   flowVGridz(_flowVGridz), flowVr(_flowVr),flowVz(_flowVz), flowVt(_flowVt),
   nR_Dens(_nR_Dens), nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),ni(_ni),
           nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr), TempGridz(_TempGridz),
           ti(_ti),te(_te),background_Z(_background_Z), background_amu(_background_amu),
   nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridR(_BfieldGridR), 
    BfieldGridZ(_BfieldGridZ),BfieldR(_BfieldR), BfieldZ(_BfieldZ), BfieldT(_BfieldT) {} 

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 

	    if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0)
        { 
        float pi = 3.14159265;    
		float nu_friction = 0.0;
		float nu_deflection = 0.0;
		float nu_parallel = 0.0;
		float nu_energy = 0.0;
		float flowVelocity[3]= {0.0, 0.0, 0.0};
		float relativeVelocity[3] = {0.0, 0.0, 0.0};
		float velocityCollisions[3];	
		float velocityRelativeNorm;	
		float parallel_direction[3];
		float perp_direction1[3];
		float perp_direction2[3];
		float parallel_contribution;
		float dv_perp1[3] = {0.0};
		float dv_perp2[3] = {0.0};
#if USE3DTETGEOM
        interp3dVector (&flowVelocity[0], particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_flowV,nY_flowV,nZ_flowV,
                flowVGridr,flowVGridy,flowVGridz,flowVr,flowVz,flowVt);
#else    
                interp2dVector(&flowVelocity[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_flowV,nZ_flowV,
                                       flowVGridr,flowVGridz,flowVr,flowVz,flowVt);
#endif
            relativeVelocity[0] = particlesPointer->vx[indx] - flowVelocity[0];
        	relativeVelocity[1] = particlesPointer->vy[indx] - flowVelocity[1];
        	relativeVelocity[2] = particlesPointer->vz[indx] - flowVelocity[2];
        	velocityRelativeNorm = vectorNorm(relativeVelocity);
            //std::cout << "particle velocity " << particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx] << std::endl;
            //std::cout << "flow velocity " << flowVelocity[0] << " " << flowVelocity[1] << " " << flowVelocity[2] << std::endl;
            //std::cout << "relative velocity and norm " << relativeVelocity[0] << " " << relativeVelocity[1] << " " << relativeVelocity[2]<< " " << velocityRelativeNorm << std::endl;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
        //int plus_minus1 = 1;//floor(curand_uniform(&particlesPointer->streams_collision1[indx]) + 0.5)*2 -1;
		//int plus_minus2 = 1;//floor(curand_uniform(&particlesPointer->streams_collision2[indx]) + 0.5)*2 -1;
		//int plus_minus3 = 1;//floor(curand_uniform(&particlesPointer->streams_collision3[indx]) + 0.5)*2 -1;
            float n1 = curand_normal(&state[indx]);
            float n2 = curand_normal(&state[indx]);
            float xsi = curand_uniform(&state[indx]);
            particlesPointer->test[indx] = xsi;
#else
	    //std::uniform_real_distribution<float> dist(0.0, 1.0);
        //int plus_minus1 = floor(dist(state[indx]) + 0.5)*2 - 1;
		//int plus_minus2 = floor(dist(state[indx]) + 0.5)*2 - 1;
		//int plus_minus3 = floor(dist(state[indx]) + 0.5)*2 - 1;
            std::normal_distribution<double> distribution(0.0,1.0);
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            float n2 = distribution(state[indx]);
            float xsi = dist(state[indx]);
#endif
#else
#if __CUDACC__
            //float plus_minus1 = floor(curand_uniform(&state[3]) + 0.5)*2-1;
            //float plus_minus2 = floor(curand_uniform(&state[4]) + 0.5)*2-1;
            //float plus_minus3 = floor(curand_uniform(&state[5]) + 0.5)*2-1;
            float n2 = curand_normal(&state[indx]);
            float xsi = curand_uniform(&state[indx]);
#else
            std::normal_distribution<double> distribution(0.0,1.0);
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            float n2 = distribution(state[indx]);
            float xsi = dist(state[indx]);
            //float plus_minus1 = floor(dist(state[3]) + 0.5)*2 - 1;
            //float plus_minus2 = floor(dist(state[4]) + 0.5)*2 - 1;
            //float plus_minus3 = floor(dist(state[5]) + 0.5)*2 - 1;
#endif
#endif 
        
        //std::cout << "random numbers n2 and xsi " << n2 << " " << xsi  << std::endl;
        
//        std::cout << "particle v " << particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx] << std::endl;
//        std::cout << "speed " << velocityRelativeNorm << std::endl;
        

		getSlowDownFrequencies (nu_friction,nu_deflection,nu_parallel, nu_energy,
        particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
        particlesPointer->vx[indx],particlesPointer->vy[indx],particlesPointer->vz[indx],
        particlesPointer->charge[indx],particlesPointer->amu[indx],
                  nR_flowV,  nZ_flowV, flowVGridr,
        flowVGridz,flowVr,
        flowVz,flowVt,
         nR_Dens, nZ_Dens,DensGridr,
        DensGridz,ni, nR_Temp,  nZ_Temp,
        TempGridr, TempGridz,ti,te, background_Z, background_amu);

		getSlowDownDirections(parallel_direction, perp_direction1,perp_direction2,
        particlesPointer->x[indx],particlesPointer->y[indx],particlesPointer->z[indx],
        particlesPointer->vx[indx],particlesPointer->vy[indx],particlesPointer->vz[indx],
                  nR_flowV,nY_flowV,  nZ_flowV, flowVGridr,flowVGridy,
        flowVGridz,flowVr,
        flowVz,flowVt,
                
    nR_Bfield,
    nZ_Bfield,
    BfieldGridR,
    BfieldGridZ,
    BfieldR,
    BfieldZ,
    BfieldT
                );
		//parallel_contribution = (1.0-nu_friction*dt + plus_minus1*sqrt(nu_parallel*dt));
		//dv_perp1[0] = perp_direction1[0]*plus_minus2*sqrt(nu_deflection/2*dt);
		//dv_perp1[1] = perp_direction1[1]*plus_minus2*sqrt(nu_deflection/2*dt);
		//dv_perp1[2] = perp_direction1[2]*plus_minus2*sqrt(nu_deflection/2*dt);
        //        dv_perp2[0] = perp_direction2[0]*plus_minus3*sqrt(nu_deflection/2*dt);
        //        dv_perp2[1] = perp_direction2[1]*plus_minus3*sqrt(nu_deflection/2*dt);
        //        dv_perp2[2] = perp_direction2[2]*plus_minus3*sqrt(nu_deflection/2*dt);
        //std::cout << "parallel direction " << parallel_direction[0] << " " << parallel_direction[1] << " " <<parallel_direction[2] << std::endl;
        //std::cout << "perp direction1 " << perp_direction1[0] << " " << perp_direction1[1] << " " <<perp_direction1[2] << std::endl;
        //std::cout << "perp direction2 " << perp_direction2[0] << " " << perp_direction2[1] << " " <<perp_direction2[2] << std::endl;
        float coeff_par = 1.0 + n1*sqrt(nu_parallel*dt);
        float coeff_perp =  abs(n2)*sqrt(nu_deflection*dt*0.5);
        float coeff_Energy = 1.0-nu_energy*dt*0.5;
        float cosXsi = cos(2.0*pi*xsi);
        float sinXsi = sin(2.0*pi*xsi);
		velocityCollisions[0] = coeff_par*parallel_direction[0] + coeff_perp*(perp_direction1[0]*cosXsi + perp_direction2[0]*sinXsi);
		velocityCollisions[1] = coeff_par*parallel_direction[1] + coeff_perp*(perp_direction1[1]*cosXsi + perp_direction2[1]*sinXsi);
		velocityCollisions[2] = coeff_par*parallel_direction[2] + coeff_perp*(perp_direction1[2]*cosXsi + perp_direction2[2]*sinXsi);
float velocityCollisionsNorm = vectorNorm(velocityCollisions);
		//velocityCollisions[0] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[0]*parallel_contribution
 		//			+ dv_perp1[0] + dv_perp2[0]);
		//velocityCollisions[1] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[1]*parallel_contribution                                         + dv_perp1[1] + dv_perp2[1]);
		//velocityCollisions[2] = velocityRelativeNorm*(1-nu_energy*dt)*(parallel_direction[2]*parallel_contribution                                         + dv_perp1[2] + dv_perp2[2]);
/*
        if (particlesPointer->charge[indx] > 0)
{
        std::cout << "nu friction, parallel perp energy " << nu_friction << " " << nu_parallel << " " <<nu_deflection << " " << nu_energy << std::endl;
        std::cout << "parallel direction " << parallel_direction[0] << " " << parallel_direction[1] << " " <<parallel_direction[2] << std::endl;
        std::cout << "parallel contribution " << parallel_contribution << std::endl;

    std::cout << "relative speed then dv "<< velocityRelativeNorm << " " << velocityCollisions[0] + flowVelocity[0] - particlesPointer->vx[indx] << " " << velocityCollisions[1] + flowVelocity[1] - particlesPointer->vy[indx]  << " " <<velocityCollisions[2] + flowVelocity[2] - particlesPointer->vz[indx]  << std::endl;
}
*/
/*        
    std::cout << "particle velocity "<< particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx]    << std::endl;
    std::cout << "flow velocity "<< flowVelocity[0] << " " << flowVelocity[1] << " " << flowVelocity[2]    << std::endl;
    std::cout << "relative velocity "<< relativeVelocity[0] << " " << relativeVelocity[1] << " " << relativeVelocity[2]    << std::endl;
    std::cout << "dv collisions "<< velocityCollisions[0] << " " << velocityCollisions[1] << " " << velocityCollisions[2]    << std::endl;
    std::cout << "relative speed then dv "<< velocityRelativeNorm << " " << velocityCollisions[0] - flowVelocity[0] - particlesPointer->vx[indx] << " " << velocityCollisions[1] - flowVelocity[1] - particlesPointer->vy[indx]  << " " <<velocityCollisions[2] + flowVelocity[2] - particlesPointer->vz[indx]  << std::endl;
*/   
    //particlesPointer->vx[indx] = velocityCollisions[0] + flowVelocity[0]; 
	//	particlesPointer->vy[indx] = velocityCollisions[1] + flowVelocity[1];
	//	particlesPointer->vz[indx] = velocityCollisions[2] + flowVelocity[2];   	
        particlesPointer->vx[indx] = coeff_Energy*velocityRelativeNorm*velocityCollisions[0]/velocityCollisionsNorm; 
		particlesPointer->vy[indx] = coeff_Energy*velocityRelativeNorm*velocityCollisions[1]/velocityCollisionsNorm;
		particlesPointer->vz[indx] = coeff_Energy*velocityRelativeNorm*velocityCollisions[2]/velocityCollisionsNorm;   	
    //std::cout << "particle velocity "<< particlesPointer->vx[indx] << " " << particlesPointer->vy[indx] << " " << particlesPointer->vz[indx]    << std::endl;
    
	}
    	}
     
};

#endif
