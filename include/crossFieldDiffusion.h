#ifndef _CFDIFFUSION_
#define _CFDIFFUSION_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

/* How do particles move perpendicular to the B field */
struct crossFieldDiffusion { 
  /* control flow */
    Flags* flags;
    /* particles we are operating on - changing their positions */
    Particles *particlesPointer;
    /* prior to time adaptivity dt was a constant, now it can adapt */
    const gitr_precision dt;
    /* how large is the step */
    /* most unknown thing is plasma physics! right here */
	const gitr_precision diffusionCoefficient;
    int nR_Bfield;
    int nZ_Bfield;
    /* 2D grid - describes Bfield spacial domain */
    gitr_precision * BfieldGridRDevicePointer;
    gitr_precision * BfieldGridZDevicePointer;
    /* particles are diffusing perpendicular to */
    /* ( R0, Z0, T0 ) ----> particle 0? No! map each particle into where it belongs spacially
       in that vector field  */
    gitr_precision * BfieldRDevicePointer;
    gitr_precision * BfieldZDevicePointer;
    gitr_precision * BfieldTDevicePointer;
#if __CUDACC__
        curandState *state;
#else
        std::mt19937 *state;
#endif
    crossFieldDiffusion(Flags* _flags, Particles *_particlesPointer, gitr_precision _dt,
#if __CUDACC__
                            curandState *_state,
#else
                                            std::mt19937 *_state,
#endif
            gitr_precision _diffusionCoefficient,
            int _nR_Bfield, int _nZ_Bfield,
            gitr_precision * _BfieldGridRDevicePointer,gitr_precision * _BfieldGridZDevicePointer,
            gitr_precision * _BfieldRDevicePointer,gitr_precision * _BfieldZDevicePointer,
            gitr_precision * _BfieldTDevicePointer)
      : flags(_flags), particlesPointer(_particlesPointer),
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

/* Monte Carlo solution to diffusion equation - we need this tested */
/* semi-non-deterministic test - tolerance type test */
/* straight field lines */
/* diffusion equation - start particles at a known point. */
CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const { 

	    if(particlesPointer->hitWall[indx] == 0.0)
        {
           if(particlesPointer->charge[indx] > 0.0)
           { 
       
	        gitr_precision perpVector[3]= {0, 0, 0};
	        gitr_precision B[3] = {0.0,0.0,0.0};
            gitr_precision Bmag = 0.0;
		gitr_precision B_unit[3] = {0.0, 0.0, 0.0};
		gitr_precision phi_random;
		gitr_precision norm;
		gitr_precision step;
		gitr_precision dt_step = dt;

                if (flags->USE_ADAPTIVE_DT) {
	          if(particlesPointer->advance[indx])
		  {
                    dt_step = particlesPointer->dt[indx];
		  }
		  else
		  {
			  dt_step = 0.0;
		  }
		}
    gitr_precision x0 = particlesPointer->xprevious[indx];
    gitr_precision y0 = particlesPointer->yprevious[indx];
    gitr_precision z0 = particlesPointer->zprevious[indx];
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
        	gitr_precision r3 = curand_uniform(&state[indx]);
#else
        	std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
        	gitr_precision r3=dist(state[indx]);
        	gitr_precision r4=dist(state[indx]);
#endif 
#else
#if __CUDACC__
            gitr_precision r3 = curand_uniform(&state[2]);
#else
            std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
            gitr_precision r3=dist(state[2]);
#endif
#endif
            /* magnitude of spacial step for 1 particle? */
            /* m^2 / sec units for diffusionCoefficient */
		step = std::sqrt(4*diffusionCoefficient*dt_step);
#if USEPERPDIFFUSION > 1
    gitr_precision plus_minus1 = floor(r4 + 0.5)*2 - 1;
    gitr_precision h = 0.001;
//    gitr_precision k1x = B_unit[0]*h;
//    gitr_precision k1y = B_unit[1]*h;
//    gitr_precision k1z = B_unit[2]*h;
    gitr_precision x_plus = x0+B_unit[0]*h;
    gitr_precision y_plus = y0+B_unit[1]*h;
    gitr_precision z_plus = z0+B_unit[2]*h;
//    gitr_precision x_minus = x0-B_unit[0]*h;
//    gitr_precision y_minus = y0-B_unit[1]*h;
//    gitr_precision z_minus = z0-B_unit[2]*h;
//     x_plus = x0+k1x;
//     y_plus = y0+k1y;
//     z_plus = z0+k1z;
//     x_minus = x0-k1x;
//     y_minus = y0-k1y;
//     z_minus = z0-k1z;
    //std::cout << "pos plus " << x_plus << " " << y_plus << " " << z_plus << std::endl;
    //std::cout << "pos minus " << x_minus << " " << y_minus << " " << z_minus << std::endl;
    gitr_precision B_plus[3] = {0.0f};
        interp2dVector(&B_plus[0],x_plus,y_plus,z_plus,nR_Bfield,nZ_Bfield,
                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
        gitr_precision Bmag_plus = std::sqrt(B_plus[0]*B_plus[0] + B_plus[1]*B_plus[1] + B_plus[2]*B_plus[2]);
//    gitr_precision k2x = B_plus[0]*h;
//    gitr_precision k2y = B_plus[1]*h;
//    gitr_precision k2z = B_plus[2]*h;
//   gitr_precision xNew = x0+0.5*(k1x+k2x); 
//   gitr_precision yNew = y0+0.5*(k1y+k2y); 
//   gitr_precision zNew = z0+0.5*(k1z+k2z); 
//   std::cout <<"pps new plus " << xNew << " " << yNew << " " << zNew << std::endl;
//    gitr_precision B_minus[3] = {0.0f};
//        interp2dVector(&B_minus[0],x_minus,y_minus,z_minus,nR_Bfield,nZ_Bfield,
//                               BfieldGridRDevicePointer,BfieldGridZDevicePointer,BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer);
//        gitr_precision Bmag_minus = std::sqrt(B_minus[0]*B_minus[0] + B_minus[1]*B_minus[1] + B_minus[2]*B_minus[2]);
//    gitr_precision k2x_minus = -B_minus[0]*h/Bmag_minus;
//    gitr_precision k2y_minus = -B_minus[1]*h/Bmag_minus;
//    gitr_precision k2z_minus = -B_minus[2]*h/Bmag_minus;
//   gitr_precision xNew_minus = x0+0.5*(k1x+k2x); 
//   gitr_precision yNew_minus = y0+0.5*(k1y+k2y); 
//   gitr_precision zNew_minus = z0+0.5*(k1z+k2z); 
//   std::cout <<"pps new minus " << xNew_minus << " " << yNew_minus << " " << zNew_minus << std::endl;
    
    gitr_precision B_deriv1[3] = {0.0f};
//    gitr_precision B_deriv2[3] = {0.0f};
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
    //gitr_precision pos_deriv1[3] = {0.0f};
    //gitr_precision pos_deriv2[3] = {0.0f};
    //pos_deriv1[0] = (xNew-x0)/(h);
    //pos_deriv1[1] = (yNew-y0)/(h);
    //pos_deriv1[2] = (zNew-z0)/(h);
    //pos_deriv2[0] = (xNew - 2*x0 + xNew_minus)/(h*h);
    //pos_deriv2[1] = (yNew - 2*y0 + yNew_minus)/(h*h);
    //pos_deriv2[2] = (zNew - 2*z0 + zNew_minus)/(h*h);
    //std::cout << "pos_deriv1 " << pos_deriv1[0] << " " <<  pos_deriv1[1]<< " " <<  pos_deriv1[2]<< " " <<std::endl;
    //gitr_precision deriv_cross[3] = {0.0};
    //vectorCrossProduct(B_deriv1, B_deriv2, deriv_cross);
    gitr_precision denom = vectorNorm(B_deriv1);
    //gitr_precision norm_cross = vectorNorm(deriv_cross);
    //std::cout << "deriv_cross " << deriv_cross[0] << " " <<  deriv_cross[1]<< " " <<  deriv_cross[2]<< " " <<std::endl;
    //std::cout << "denome and norm_cross " << denom << " " << norm_cross << std::endl;
    gitr_precision R = 1.0e4;
    if((std::abs(denom) > 1e-10) & (std::abs(denom) < 1e10) )
    {
      R = Bmag/denom;
    }
    //std::cout << "Radius of curvature"<< R <<std::endl;
    gitr_precision initial_guess_theta = 3.14159265359*0.5;
    gitr_precision eps = 0.01;
    gitr_precision error = 2.0;
    gitr_precision s = step;
    gitr_precision drand = r3;
    gitr_precision theta0 = initial_guess_theta;
    gitr_precision theta1 = 0.0;
    gitr_precision f = 0.0;
    gitr_precision f_prime = 0.0;
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
    gitr_precision y_dir[3] = {0.0};
    vectorCrossProduct(B, B_deriv1, y_dir);
    gitr_precision x_comp = s*std::cos(theta0);
    gitr_precision y_comp = s*std::sin(theta0);
    gitr_precision x_transform = x_comp*perpVector[0] + y_comp*y_dir[0];
    gitr_precision y_transform = x_comp*perpVector[1] + y_comp*y_dir[1];
    gitr_precision z_transform = x_comp*perpVector[2] + y_comp*y_dir[2];
    //gitr_precision R = 1.2;
    //gitr_precision r_in = R-step;
    //gitr_precision r_out = R-step;
    //gitr_precision A_in = R*R-r_in*r_in;
    //gitr_precision A_out = r_out*r_out-R*R;
    //gitr_precision theta[100] = {0.0f};
    //gitr_precision f_theta[100] = {0.0f};
    //gitr_precision cdf[100] = {0.0f};
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
    //gitr_precision minTheta = 0;
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
