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
      : particlesPointer(_particlesPointer),
        dt(_dt),
        diffusionCoefficient(_diffusionCoefficient),
        nR_Bfield(_nR_Bfield),
        nZ_Bfield(_nZ_Bfield),
        BfieldGridRDevicePointer(_BfieldGridRDevicePointer),
        BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
        BfieldRDevicePointer(_BfieldRDevicePointer),
        BfieldZDevicePointer(_BfieldZDevicePointer),
        BfieldTDevicePointer(_BfieldTDevicePointer),
        state(_state) {
  }

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
    float x0 = particlesPointer->xprevious[indx];
    float y0 = particlesPointer->yprevious[indx];
    float z0 = particlesPointer->zprevious[indx];
    //std::cout << "initial position " << x0 << " " << y0 << " " << z0 << std::endl;
        interp2dVector(&B[0],particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_Bfield,nZ_Bfield,
                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
        Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        B_unit[0] = B[0]/Bmag;
        B_unit[1] = B[1]/Bmag;
        B_unit[2] = B[2]/Bmag;
    //std::cout << "B " << B[0] << " " <<  B[1]<< " " <<  B[2]<< " " <<std::endl;
    //std::cout << "B_unit " << B_unit[0] << " " <<  B_unit[1]<< " " <<  B_unit[2]<< " " <<std::endl;
#if PARTICLESEEDS > 0
#ifdef __CUDACC__
        	float r3 = curand_uniform(&state[indx]);
#else
        	std::uniform_real_distribution<float> dist(0.0, 1.0);
        	float r3=dist(state[indx]);
        	float r4=dist(state[indx]);
#endif 
#else
#if __CUDACC__
            float r3 = curand_uniform(&state[2]);
#else
            std::uniform_real_distribution<float> dist(0.0, 1.0);
            float r3=dist(state[2]);
#endif
#endif
		step = std::sqrt(6*diffusionCoefficient*dt);
#if USEPERPDIFFUSION > 1
    float plus_minus1 = floor(r4 + 0.5)*2 - 1;
    float h = 0.001;
//    float k1x = B_unit[0]*h;
//    float k1y = B_unit[1]*h;
//    float k1z = B_unit[2]*h;
    float x_plus = x0+B_unit[0]*h;
    float y_plus = y0+B_unit[1]*h;
    float z_plus = z0+B_unit[2]*h;
//    float x_minus = x0-B_unit[0]*h;
//    float y_minus = y0-B_unit[1]*h;
//    float z_minus = z0-B_unit[2]*h;
//     x_plus = x0+k1x;
//     y_plus = y0+k1y;
//     z_plus = z0+k1z;
//     x_minus = x0-k1x;
//     y_minus = y0-k1y;
//     z_minus = z0-k1z;
    //std::cout << "pos plus " << x_plus << " " << y_plus << " " << z_plus << std::endl;
    //std::cout << "pos minus " << x_minus << " " << y_minus << " " << z_minus << std::endl;
    float B_plus[3] = {0.0f};
        interp2dVector(&B_plus[0],x_plus,y_plus,z_plus,nR_Bfield,nZ_Bfield,
                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
        float Bmag_plus = std::sqrt(B_plus[0]*B_plus[0] + B_plus[1]*B_plus[1] + B_plus[2]*B_plus[2]);
//    float k2x = B_plus[0]*h;
//    float k2y = B_plus[1]*h;
//    float k2z = B_plus[2]*h;
//   float xNew = x0+0.5*(k1x+k2x); 
//   float yNew = y0+0.5*(k1y+k2y); 
//   float zNew = z0+0.5*(k1z+k2z); 
//   std::cout <<"pps new plus " << xNew << " " << yNew << " " << zNew << std::endl;
//    float B_minus[3] = {0.0f};
//        interp2dVector(&B_minus[0],x_minus,y_minus,z_minus,nR_Bfield,nZ_Bfield,
//                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
//        float Bmag_minus = std::sqrt(B_minus[0]*B_minus[0] + B_minus[1]*B_minus[1] + B_minus[2]*B_minus[2]);
//    float k2x_minus = -B_minus[0]*h/Bmag_minus;
//    float k2y_minus = -B_minus[1]*h/Bmag_minus;
//    float k2z_minus = -B_minus[2]*h/Bmag_minus;
//   float xNew_minus = x0+0.5*(k1x+k2x); 
//   float yNew_minus = y0+0.5*(k1y+k2y); 
//   float zNew_minus = z0+0.5*(k1z+k2z); 
//   std::cout <<"pps new minus " << xNew_minus << " " << yNew_minus << " " << zNew_minus << std::endl;
    
    float B_deriv1[3] = {0.0f};
//    float B_deriv2[3] = {0.0f};
    //std::cout << "B_plus " << B_plus[0] << " " <<  B_plus[1]<< " " <<  B_plus[2]<< " " <<std::endl;
    //std::cout << "B_minus " << B_minus[0] << " " <<  B_minus[1]<< " " <<  B_minus[2]<< " " <<std::endl;
    B_deriv1[0] = (B_plus[0] - B[0])/(h);
    B_deriv1[1] = (B_plus[1] - B[1])/(h);
    B_deriv1[2] = (B_plus[2] - B[2])/(h);
    
   // B_deriv1[0] = (B_plus[0] - B_minus[0])/(2*h);
   // B_deriv1[1] = (B_plus[1] - B_minus[1])/(2*h);
   // B_deriv1[2] = (B_plus[2] - B_minus[2])/(2*h);
    //std::cout << "B_deriv1 " << B_deriv1[0] << " " <<  B_deriv1[1]<< " " <<  B_deriv1[2]<< " " <<std::endl;
    //std::cout << "Bderiv2 " << B_deriv2[0] << " " <<  B_deriv2[1]<< " " <<  B_deriv2[2]<< " " <<std::endl;
    //float pos_deriv1[3] = {0.0f};
    //float pos_deriv2[3] = {0.0f};
    //pos_deriv1[0] = (xNew-x0)/(h);
    //pos_deriv1[1] = (yNew-y0)/(h);
    //pos_deriv1[2] = (zNew-z0)/(h);
    //pos_deriv2[0] = (xNew - 2*x0 + xNew_minus)/(h*h);
    //pos_deriv2[1] = (yNew - 2*y0 + yNew_minus)/(h*h);
    //pos_deriv2[2] = (zNew - 2*z0 + zNew_minus)/(h*h);
    //std::cout << "pos_deriv1 " << pos_deriv1[0] << " " <<  pos_deriv1[1]<< " " <<  pos_deriv1[2]<< " " <<std::endl;
    //float deriv_cross[3] = {0.0};
    //vectorCrossProduct(B_deriv1, B_deriv2, deriv_cross);
    float denom = vectorNorm(B_deriv1);
    //float norm_cross = vectorNorm(deriv_cross);
    //std::cout << "deriv_cross " << deriv_cross[0] << " " <<  deriv_cross[1]<< " " <<  deriv_cross[2]<< " " <<std::endl;
    //std::cout << "denome and norm_cross " << denom << " " << norm_cross << std::endl;
    float R = 1.0e4;
    if((std::abs(denom) > 1e-10) & (std::abs(denom) < 1e10) )
    {
      R = Bmag/denom;
    }
    //std::cout << "Radius of curvature"<< R <<std::endl;
    float initial_guess_theta = 3.14159265359*0.5;
    float eps = 0.01;
    float error = 2.0;
    float s = step;
    float drand = r3;
    float theta0 = initial_guess_theta;
    float theta1 = 0.0;
    float f = 0.0;
    float f_prime = 0.0;
    int nloops = 0;
    if(R > 1.0e-4)
    {
    while ((error > eps)&(nloops<10))
    {
        f = (2*R*theta0-s*sin(theta0))/(2*3.14159265359*R) - drand;
        f_prime = (2*R-s*cos(theta0))/(2*3.14159265359*R); 
        theta1 = theta0 - f/f_prime;
         error = abs(theta1-theta0);
         theta0=theta1;
         nloops++;
         //std::cout << " R rand and theta "<<R << " " <<  drand << " " << theta0 << std::endl;
    }
    if(nloops > 9)
    {
      theta0 = 2*3.14159265359*drand;

    }
    }
    else
    {R = 1.0e-4;
      theta0 = 2*3.14159265359*drand;
    }
         //std::cout << "out of newton"<< std::endl;
      

    if(plus_minus1 < 0)
    {
      theta0 = 2*3.14159265359-theta0; 
    }
perpVector[0] = B_deriv1[0];
perpVector[1] = B_deriv1[1];
perpVector[2] = B_deriv1[2];
		norm = std::sqrt(perpVector[0]*perpVector[0] + perpVector[1]*perpVector[1] + perpVector[2]*perpVector[2]);
		perpVector[0] = perpVector[0]/norm;
		perpVector[1] = perpVector[1]/norm;
		perpVector[2] = perpVector[2]/norm;
    float y_dir[3] = {0.0};
    vectorCrossProduct(B, B_deriv1, y_dir);
    float x_comp = s*std::cos(theta0);
    float y_comp = s*std::sin(theta0);
    float x_transform = x_comp*perpVector[0] + y_comp*y_dir[0];
    float y_transform = x_comp*perpVector[1] + y_comp*y_dir[1];
    float z_transform = x_comp*perpVector[2] + y_comp*y_dir[2];
    //float R = 1.2;
    //float r_in = R-step;
    //float r_out = R-step;
    //float A_in = R*R-r_in*r_in;
    //float A_out = r_out*r_out-R*R;
    //float theta[100] = {0.0f};
    //float f_theta[100] = {0.0f};
    //float cdf[100] = {0.0f};
    //for(int ii=0;ii<100;ii++)
    //{ theta[ii] = 3.1415/100.0*ii;
    //  f_theta[ii] = 
    //    //2*R*theta[ii] -step*std::sin(theta[ii]);}
    //    //A_in + ii/100*(A_out-A_in);}
    //    (2.0*R*step-step*step*std::cos(theta[ii]));}
    //cdf[0] = f_theta[0];
    //for(int ii=1;ii<100;ii++)
    //{cdf[ii] = cdf[ii-1]+f_theta[ii];}
    //for(int ii=0;ii<100;ii++)
    //{cdf[ii] = cdf[ii]+cdf[99];}
    //float minTheta = 0;
    //int found=0;
    //for(int ii=0;ii<100;ii++)
    //{if(r3>cdf[ii] && found < 1)
    //  {
    //    minTheta = theta[ii];
    //    found = 2;
    //  }
    //}

if(std::abs(denom) < 1.0e-8)
{
#endif
    perpVector[0] = 0.0;
    perpVector[1] = 0.0;
    perpVector[2] = 0.0;
		phi_random = 2*3.14159265*r3;
		perpVector[0] = std::cos(phi_random);
		perpVector[1] = std::sin(phi_random);
		perpVector[2] = (-perpVector[0]*B_unit[0] - perpVector[1]*B_unit[1])/B_unit[2];
                //std::cout << "perp Vector " << perpVector[0] << " " << perpVector[1] << " " << perpVector[2]<<std::endl;
		if (B_unit[2] == 0){
			perpVector[2] = perpVector[1];
			perpVector[1] = (-perpVector[0]*B_unit[0] - perpVector[2]*B_unit[2])/B_unit[1];
		}
               // std::cout << "perp Vector " << perpVector[0] << " " << perpVector[1] << " " << perpVector[2]<<std::endl;
		
		if ((B_unit[0] == 1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0) || (B_unit[0] == -1.0 && B_unit[1] ==0.0 && B_unit[2] ==0.0))
		{
			perpVector[2] = perpVector[0];
			perpVector[0] = 0;
			perpVector[1] = sin(phi_random);
                //std::cout << "perp Vector " << perpVector[0] << " " << perpVector[1] << " " << perpVector[2]<<std::endl;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==1.0 && B_unit[2] ==0.0) || (B_unit[0] == 0.0 && B_unit[1] ==-1.0 && B_unit[2] ==0.0))
		{
			perpVector[1] = 0.0;
		}
		else if ((B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==1.0) || (B_unit[0] == 0.0 && B_unit[1] ==0.0 && B_unit[2] ==-1.0))
		{
			perpVector[2] = 0;
		}
		
		norm = std::sqrt(perpVector[0]*perpVector[0] + perpVector[1]*perpVector[1] + perpVector[2]*perpVector[2]);
		perpVector[0] = perpVector[0]/norm;
		perpVector[1] = perpVector[1]/norm;
		perpVector[2] = perpVector[2]/norm;
//                //std::cout << "perp Vector " << perpVector[0] << " " << perpVector[1] << " " << perpVector[2]<<std::endl;
//		
//		step = std::sqrt(6*diffusionCoefficient*dt);
    //std::cout << "y_dir " << y_dir[0] << " " <<  y_dir[1]<< " " <<  y_dir[2]<< " " <<std::endl;
    //std::cout << "y_dir " << y_dir[0] << " " <<  y_dir[1]<< " " <<  y_dir[2]<< " " <<std::endl;
    //std::cout << "transforms " << x_transform << " " << y_transform << " " << z_transform << std::endl; 
		particlesPointer->x[indx] = particlesPointer->xprevious[indx] + step*perpVector[0];
		particlesPointer->y[indx] = particlesPointer->yprevious[indx] + step*perpVector[1];
		particlesPointer->z[indx] = particlesPointer->zprevious[indx] + step*perpVector[2];
#if USEPERPDIFFUSION > 1    
}
else
{
    particlesPointer->x[indx] = x0 + x_transform;
		particlesPointer->y[indx] = y0 + y_transform;
		particlesPointer->z[indx] = z0 + z_transform;
}
#endif
    	}
    } }
};

#endif
