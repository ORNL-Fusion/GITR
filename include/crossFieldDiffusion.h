#ifndef _CFDIFFUSION_
#define _CFDIFFUSION_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "boris.h"
#include "Particles.h"
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

CUDA_CALLABLE_MEMBER_DEVICE    
void legacy_code_block_0( Particles *particlesPointer, 
                          int indx,
                          gitr_precision *B_unit,
                          gitr_precision step,
                          gitr_precision r3 )
{
      gitr_precision perpVector[3]= {0, 0, 0};

      gitr_precision phi_random = 2*3.14159265*r3;
      perpVector[0] = std::cos(phi_random);
      perpVector[1] = std::sin(phi_random);
      perpVector[2] = 0.0;
    
      gitr_precision ez1 = B_unit[0];
      gitr_precision ez2 = B_unit[1];
      gitr_precision ez3 = B_unit[2];

      // Get perpendicular velocity unit vectors
      // this comes from a cross product of
      // (ez1,ez2,ez3)x(0,0,1)
      gitr_precision ex1 = ez2;
      gitr_precision ex2 = -ez1;
      gitr_precision ex3 = 0.0;
    
      // The above cross product will be zero for particles
      // with a pure z-directed (ez3) velocity
      // here we find those particles and get the perpendicular 
      // unit vectors by taking the cross product
      // (ez1,ez2,ez3)x(0,1,0) instead
      gitr_precision exnorm = std::sqrt(ex1*ex1 + ex2*ex2);
      if(std::abs(exnorm) < 1.0e-12)
      {
        ex1 = -ez3;
        ex2 = 0.0;
        ex3 = ez1;
      }

      // Ensure all the perpendicular direction vectors
      // ex are unit
      exnorm = std::sqrt(ex1*ex1+ex2*ex2 + ex3*ex3);
      ex1 = ex1/exnorm;
      ex2 = ex2/exnorm;
      ex3 = ex3/exnorm;
      
      // Find the second perpendicular direction 
      // by taking the cross product
      // (ez1,ez2,ez3)x(ex1,ex2,ex3)
      gitr_precision ey1 = ez2*ex3 - ez3*ex2;
      gitr_precision ey2 = ez3*ex1 - ez1*ex3;
      gitr_precision ey3 = ez1*ex2 - ez2*ex1;

      gitr_precision tmp[3] = {0.0};
      tmp[0] = ex1*perpVector[0] + ey1*perpVector[1] + ez1*perpVector[2];
      tmp[1] = ex2*perpVector[0] + ey2*perpVector[1] + ez2*perpVector[2];
      tmp[2] = ex3*perpVector[0] + ey3*perpVector[1] + ez3*perpVector[2];

      perpVector[0] = tmp[0];
      perpVector[1] = tmp[1];
      perpVector[2] = tmp[2];

      gitr_precision norm = 
      std::sqrt( perpVector[0]*perpVector[0] + 
                  perpVector[1]*perpVector[1] +
                  perpVector[2]*perpVector[2] );

      perpVector[0] = perpVector[0]/norm;
      perpVector[1] = perpVector[1]/norm;
      perpVector[2] = perpVector[2]/norm;
      particlesPointer->x[indx] = particlesPointer->xprevious[indx] + step*perpVector[0];
      particlesPointer->y[indx] = particlesPointer->yprevious[indx] + step*perpVector[1];
      particlesPointer->z[indx] = particlesPointer->zprevious[indx] + step*perpVector[2];
}

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
        int perp_diffusion;

        int cylsymm;

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
        gitr_precision * _BfieldTDevicePointer, int perp_diffusion_, int cylsymm_ )
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
        state(_state),
        perp_diffusion( perp_diffusion_ ),
        cylsymm( cylsymm_ )
        { }

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
      gitr_precision norm;
      gitr_precision step;
      gitr_precision dt_step = dt;

      if ( flags->USE_ADAPTIVE_DT ) 
      {
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
        
      interp2dVector(&B[0],particlesPointer->xprevious[indx],
                          particlesPointer->yprevious[indx],
                          particlesPointer->zprevious[indx],
                          nR_Bfield,nZ_Bfield,
                          BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                          BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer,
                          cylsymm );
        
      Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
      if(Bmag < 1.0e-12) 
      {
        B[0] = 0.0;
        B[1] = 0.0;
        B[2] = 1.0;
        Bmag = 1.0;
      }
      B_unit[0] = B[0]/Bmag;
      B_unit[1] = B[1]/Bmag;
      B_unit[2] = B[2]/Bmag;

#ifdef __CUDACC__
      gitr_precision r3 = curand_uniform(&state[indx]);
#else
      std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
      gitr_precision r3=dist(state[indx]);
      gitr_precision r4=dist(state[indx]);
#endif 

      /* magnitude of spacial step for 1 particle? */
      /* m^2 / sec units for diffusionCoefficient */
      step = std::sqrt(4*diffusionCoefficient*dt_step);
      /* Captain! Is an "else" even needed here? It doesn't appear to be so */
      if( perp_diffusion <= 1 )
      {
        legacy_code_block_0( particlesPointer, indx, B_unit, step, r3 );
      }

      else
      {
      /* Notice of code change: previously, this wass floor(r4 + 0.5)*2 - 1 */
      /* r3 was replaced with variable r4 since r4 is not declared if CUDA is activated */
      gitr_precision plus_minus1 = floor(r3 + 0.5)*2 - 1;
      gitr_precision h = 0.001;
      gitr_precision x_plus = x0+B_unit[0]*h;
      gitr_precision y_plus = y0+B_unit[1]*h;
      gitr_precision z_plus = z0+B_unit[2]*h;
     
      gitr_precision B_plus[3] = {0.0};
      
      interp2dVector(&B_plus[0],x_plus,y_plus,z_plus,nR_Bfield,nZ_Bfield,
                     BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                     BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer,
                     cylsymm );

        gitr_precision Bmag_plus = std::sqrt(B_plus[0]*B_plus[0] + B_plus[1]*B_plus[1] + B_plus[2]*B_plus[2]);
    
    gitr_precision B_deriv1[3] = {0.0};
    
    B_deriv1[0] = (B_plus[0] - B[0])/(h);
    B_deriv1[1] = (B_plus[1] - B[1])/(h);
    B_deriv1[2] = (B_plus[2] - B[2])/(h);
    
    gitr_precision denom = vectorNorm(B_deriv1);
    
    gitr_precision R = 1.0e4;
    
    if((std::abs(denom) > 1e-10) & (std::abs(denom) < 1e10) )
    {
      R = Bmag/denom;
    }
    
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
    {
      R = 1.0e-4;
      theta0 = 2*3.14159265359*drand;
    }

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

    if(std::abs(denom) < 1.0e-8)
    {
//#endif
    /* Captain! Turn this into a function. Move stuff from the block into the comment then
       move it all and butcher the control flow. Define and call the lambda function here */
    //f( perpVector, particlesPointer, indx, B_unit, step )
    legacy_code_block_0( particlesPointer, indx, B_unit, step, r3 );

  /* Captain! */
    }
else
{
    particlesPointer->x[indx] = x0 + x_transform;
    particlesPointer->y[indx] = y0 + y_transform;
    particlesPointer->z[indx] = z0 + z_transform;
}
      }
//#endif
    }
    } }
};

#endif
