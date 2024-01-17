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
#include "interp2d.hpp"
#include "boris.h"
#include "array.h"
#include "flags.hpp"
#include <iomanip>

#ifdef __CUDACC__
#include <thrust/random.h>
#else
#include <random>
#endif
#include <fenv.h>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

/* drop-in replacement for interp2dCombined() */
CUDA_CALLABLE_MEMBER
gitr_precision intorp( 
  gitr_precision x, 
  gitr_precision y, gitr_precision z,
  int nx,
  int nz,
  gitr_precision* gridx,
  gitr_precision* gridz,
  gitr_precision* data,
  int cylsymm )
{
    gitr_precision fxz = 0.0;
    gitr_precision fx_z1 = 0.0;
    gitr_precision fx_z2 = 0.0; 
    if(nx*nz == 1)
    {
        fxz = data[0];
    }
    else{
    gitr_precision dim1;
     if( cylsymm )
     {
    dim1 = std::sqrt(x*x + y*y);
    }
    else
    {
    dim1 = x;
    }
    gitr_precision d_dim1 = gridx[1] - gridx[0];
    gitr_precision dz = gridz[1] - gridz[0];
    int i = std::floor((dim1 - gridx[0])/d_dim1);//addition of 0.5 finds nearest gridpoint
    int j = std::floor((z - gridz[0])/dz);
    
    //gitr_precision interp_value = data[i + j*nx];
    if (i < 0) i=0;
    if (j < 0) j=0;
    if (i >=nx-1 && j>=nz-1)
    {
        fxz = data[nx-1+(nz-1)*nx];
    }
    else if (i >=nx-1)
    {
        fx_z1 = data[nx-1+j*nx];
        fx_z2 = data[nx-1+(j+1)*nx];
        fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
    }
    else if (j >=nz-1)
    {
        fx_z1 = data[i+(nz-1)*nx];
        fx_z2 = data[i+(nz-1)*nx];
        fxz = ((gridx[i+1]-dim1)*fx_z1+(dim1 - gridx[i])*fx_z2)/d_dim1;
        
    }
    else
    {
      fx_z1 = ((gridx[i+1]-dim1)*data[i+j*nx] + (dim1 - gridx[i])*data[i+1+j*nx])/d_dim1;
      fx_z2 = ((gridx[i+1]-dim1)*data[i+(j+1)*nx] + (dim1 - gridx[i])*data[i+1+(j+1)*nx])/d_dim1; 
      fxz = ((gridz[j+1]-z)*fx_z1+(z - gridz[j])*fx_z2)/dz;
      /* break 0 */
    }
    }
  /* Captain! new code begin */
  long long unsigned int dims[ 2 ] = { nz, nx };
  double min_range_init[ 2 ] = { gridz[ 0 ], gridx[ 0 ] };
  double max_range_init[ 2 ] = { gridz[ dims[ 0 ] - 1 ], gridx[ dims[ 1 ] - 1 ] };
  int n_dims_init = 2;

  interpolated_field< double >
  ifield( data, dims, max_range_init, min_range_init, n_dims_init );

  double coordinates[ 2 ] = { z, x };

  double value = ifield( coordinates );

  double diff = value - fxz;
  std::cout << std::setprecision(10) 
            << "value: " << value << " fxz: " << fxz << " diff: " << diff << std::endl;

  /* new code end */


    return fxz;
}


CUDA_CALLABLE_MEMBER
void getSlowDownFrequencies ( gitr_precision& nu_friction, gitr_precision& nu_deflection, gitr_precision& nu_parallel,gitr_precision& nu_energy, gitr_precision x, gitr_precision y,gitr_precision z, gitr_precision vx, gitr_precision vy, gitr_precision vz,gitr_precision charge, gitr_precision amu,
    int nR_flowV,
    int nZ_flowV,
    gitr_precision* flowVGridr,
    gitr_precision* flowVGridz,
    gitr_precision* flowVr,
    gitr_precision* flowVz,
    gitr_precision* flowVt,
    int nR_Dens,
    int nZ_Dens,
    gitr_precision* DensGridr,
    gitr_precision* DensGridz,
    gitr_precision* ni,
    int nR_Temp,
    int nZ_Temp,
    gitr_precision* TempGridr,
    gitr_precision* TempGridz,
    gitr_precision* ti, gitr_precision* te,gitr_precision background_Z, gitr_precision background_amu,
    int nR_Bfield, int nZ_Bfield,
    gitr_precision* BfieldGridR ,gitr_precision* BfieldGridZ ,
    gitr_precision* BfieldR ,gitr_precision* BfieldZ ,
    gitr_precision* BfieldT,gitr_precision &T_background,
    int flowv_interp, int cylsymm, int field_aligned_values, bool use_sheath_density, gitr_precision f_psi )
{
  //int feenableexcept(FE_INVALID | FE_OVERFLOW); //enables trapping of the floating-point exceptions
  gitr_precision Q = 1.60217662e-19;
  gitr_precision EPS0 = 8.854187e-12;
  gitr_precision pi = 3.14159265;
  gitr_precision MI = 1.6737236e-27;	
  gitr_precision ME = 9.10938356e-31;
        
  gitr_precision te_eV = intorp(x,y,z,nR_Temp,nZ_Temp,TempGridr,TempGridz,te, 
                                          cylsymm );
  gitr_precision ti_eV = intorp(x,y,z,nR_Temp,nZ_Temp,TempGridr,TempGridz,ti, 
                         cylsymm );

  T_background = ti_eV;
  gitr_precision density = interp2dCombined( x,y,z,nR_Dens,nZ_Dens,DensGridr,DensGridz,ni,
                                             cylsymm );

  if( use_sheath_density )
  {
    density = density*f_psi;
  }
  gitr_precision flowVelocity[3]= {0.0};
  gitr_precision relativeVelocity[3] = {0.0, 0.0, 0.0};
  gitr_precision velocityNorm = 0.0;
  gitr_precision lam_d;
  gitr_precision lam;
  gitr_precision gam_electron_background;
  gitr_precision gam_ion_background;
  gitr_precision a_electron = 0.0;
  gitr_precision a_ion = 0.0;
  gitr_precision xx;
  gitr_precision psi_prime;
  gitr_precision psi_psiprime;
  gitr_precision psi;
  gitr_precision xx_e;
  gitr_precision psi_prime_e;
  gitr_precision psi_psiprime_e;
  gitr_precision psi_psiprime_psi2x = 0.0;
  gitr_precision psi_psiprime_psi2x_e = 0.0;
  gitr_precision psi_e;
  gitr_precision nu_0_i;
  gitr_precision nu_0_e;
  gitr_precision nu_friction_i;
  gitr_precision nu_deflection_i;
  gitr_precision nu_parallel_i;
  gitr_precision nu_energy_i;
  gitr_precision nu_friction_e;
  gitr_precision nu_deflection_e;
  gitr_precision nu_parallel_e;
  gitr_precision nu_energy_e;
                
  if( flowv_interp == 3 )
  {
    exit( 0 );
    /*
  interp3dVector (&flowVelocity[0], x,y,z,nR_flowV,nY_flowV,nZ_flowV,
                   flowVGridr,flowVGridy,flowVGridz,flowVr,flowVz,flowVt);
    */
  }
  else if( flowv_interp < 3 )
  {
  if( field_aligned_values > 0 )
  {
  interpFieldAlignedVector(&flowVelocity[0],x,y,z,
                           nR_flowV,nZ_flowV,
                           flowVGridr,flowVGridz,flowVr,
                           flowVz,flowVt,nR_Bfield,nZ_Bfield,
                           BfieldGridR,BfieldGridZ,BfieldR,
                           BfieldZ,BfieldT, cylsymm );
  }
  else
  {
  interp2dVector(&flowVelocity[0],x,y,z,
           nR_flowV,nZ_flowV,
           flowVGridr,flowVGridz,flowVr,flowVz,flowVt, cylsymm );
  }
  }
  relativeVelocity[0] = vx - flowVelocity[0];
  relativeVelocity[1] = vy - flowVelocity[1];
  relativeVelocity[2] = vz - flowVelocity[2];
  velocityNorm = std::sqrt( relativeVelocity[0]*relativeVelocity[0] + relativeVelocity[1]*relativeVelocity[1] + relativeVelocity[2]*relativeVelocity[2]);                

  lam_d = std::sqrt(EPS0*te_eV/(density*std::pow(background_Z,2)*Q));//only one q in order to convert to J
  lam = 12.0*pi*density*std::pow(lam_d,3)/charge;
  gam_electron_background = 0.238762895*std::pow(charge,2)*std::log(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)
  gam_ion_background = 0.238762895*std::pow(charge,2)*std::pow(background_Z,2)*std::log(lam)/(amu*amu);//constant = Q^4/(MI^2*4*pi*EPS0^2)

  if(gam_electron_background < 0.0) gam_electron_background=0.0;
  if(gam_ion_background < 0.0) gam_ion_background=0.0;
  a_ion = background_amu*MI/(2*ti_eV*Q);// %q is just to convert units - no z needed
  a_electron = ME/(2*te_eV*Q);// %q is just to convert units - no z needed

  xx = std::pow(velocityNorm,2)*a_ion;
  psi_prime = 2.0*std::sqrt(xx/pi)*std::exp(-xx);
  psi_psiprime = std::erf(std::sqrt(xx));
  psi = psi_psiprime - psi_prime;
  
  xx_e = std::pow(velocityNorm,2)*a_electron;
  psi_prime_e = 1.128379*std::sqrt(xx_e);
  psi_e = 0.75225278*std::pow(xx_e,1.5);
  psi_psiprime_e = psi_e+psi_prime_e;
  psi_psiprime_psi2x_e = 1.128379*std::sqrt(xx_e)*expf(-xx_e);
  
  nu_0_i = gam_electron_background*density/std::pow(velocityNorm,3);
  nu_0_e = gam_ion_background*density/std::pow(velocityNorm,3);
  
  nu_friction_i = (1+amu/background_amu)*psi*nu_0_i;
  nu_deflection_i = 2*(psi_psiprime - psi/(2*xx))*nu_0_i;
  nu_parallel_i = psi/xx*nu_0_i;
  nu_energy_i = 2*(amu/background_amu*psi - psi_prime)*nu_0_i;
  
  nu_friction_e = (1+amu/(ME/MI))*psi_e*nu_0_e;
  nu_deflection_e = 2*(psi_psiprime_psi2x_e)*nu_0_e;
  nu_parallel_e = psi_e/xx_e*nu_0_e;
  nu_energy_e = 2*(amu/(ME/MI)*psi_e - psi_prime_e)*nu_0_e;
                    
  nu_friction = nu_friction_i ;//+ nu_friction_e;
  nu_deflection = nu_deflection_i ;//+ nu_deflection_e;
  nu_parallel = nu_parallel_i;// + nu_parallel_e;
  nu_energy = nu_energy_i;// + nu_energy_e;
    
  if(te_eV <= 0.0 || ti_eV <= 0.0)
  {
    nu_friction = 0.0;
    nu_deflection = 0.0;
    nu_parallel = 0.0;
    nu_energy = 0.0;
  }
  if(density <= 0.0)
  {
    nu_friction = 0.0;
    nu_deflection = 0.0;
    nu_parallel = 0.0;
    nu_energy = 0.0;
  }
}
CUDA_CALLABLE_MEMBER
void getSlowDownDirections2 (gitr_precision parallel_direction[], gitr_precision perp_direction1[], gitr_precision perp_direction2[],
        gitr_precision vx, gitr_precision vy, gitr_precision vz)
{
  gitr_precision v = std::sqrt(vx*vx + vy*vy + vz*vz);
  
  if(v == 0.0)
  {
    v = 1.0;
    vz = 1.0;
    vx = 0.0;
    vy = 0.0;
  }
  
  gitr_precision ez1 = vx/v;
  gitr_precision ez2 = vy/v;
  gitr_precision ez3 = vz/v;
    
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
  if(std::abs(exnorm) < 1.0e-12){
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
  
  //if(isnan(ex1) || isnan(ex2) || isnan(ex3)){
  //   printf("ex nan %f %f %f v %f", ez1, ez2, ez3,v);
  //}
  // Find the second perpendicular direction 
  // by taking the cross product
  // (ez1,ez2,ez3)x(ex1,ex2,ex3)
  gitr_precision ey1 = ez2*ex3 - ez3*ex2;
  gitr_precision ey2 = ez3*ex1 - ez1*ex3;
  gitr_precision ey3 = ez1*ex2 - ez2*ex1;
  parallel_direction[0] = ez1; 
  parallel_direction[1] = ez2;
  parallel_direction[2] = ez3;
  
  perp_direction1[0] = ex1; 
  perp_direction1[1] = ex2;
  perp_direction1[2] = ex3;
  
  perp_direction2[0] = ey1; 
  perp_direction2[1] = ey2;
  perp_direction2[2] = ey3;
}

struct coulombCollisions { 

    Particles *particlesPointer;
    gitr_precision dt;
    int nR_flowV;
    int nY_flowV;
    int nZ_flowV;
    gitr_precision* flowVGridr;
    gitr_precision* flowVGridy;
    gitr_precision* flowVGridz;
    gitr_precision* flowVr;
    gitr_precision* flowVz;
    gitr_precision* flowVt;
    int nR_Dens;
    int nZ_Dens;
    gitr_precision* DensGridr;
    gitr_precision* DensGridz;
    gitr_precision* ni;
    int nR_Temp;
    int nZ_Temp;
    gitr_precision* TempGridr;
    gitr_precision* TempGridz;
    gitr_precision* ti;
    gitr_precision* te;
    gitr_precision background_Z;
    gitr_precision background_amu;
    int nR_Bfield;
    int nZ_Bfield;
    gitr_precision * BfieldGridR;
    gitr_precision * BfieldGridZ;
    gitr_precision * BfieldR;
    gitr_precision * BfieldZ;
    gitr_precision * BfieldT;
    Flags* gitr_flags;
    gitr_precision dv[3];
#if __CUDACC__
            curandState *state;
#else
            std::mt19937 *state;
#endif

    int flowv_interp;

    int cylsymm;

    int field_aligned_values;
    int coulomb_collisions;

    coulombCollisions(Particles *_particlesPointer,gitr_precision _dt, 
#if __CUDACC__
                            curandState *_state,
#else
                            std::mt19937 *_state,
#endif
            int _nR_flowV,int _nY_flowV, int _nZ_flowV,    gitr_precision* _flowVGridr,gitr_precision* _flowVGridy,
                gitr_precision* _flowVGridz,gitr_precision* _flowVr,
                        gitr_precision* _flowVz,gitr_precision* _flowVt,
                        int _nR_Dens,int _nZ_Dens,gitr_precision* _DensGridr,
                            gitr_precision* _DensGridz,gitr_precision* _ni,int _nR_Temp, int _nZ_Temp,
                        gitr_precision* _TempGridr, gitr_precision* _TempGridz,gitr_precision* _ti,gitr_precision* _te,
                        gitr_precision _background_Z, gitr_precision _background_amu,
                        int _nR_Bfield, int _nZ_Bfield,
                        gitr_precision * _BfieldGridR ,gitr_precision * _BfieldGridZ ,
                        gitr_precision * _BfieldR ,gitr_precision * _BfieldZ ,
                 gitr_precision * _BfieldT, Flags* _gitr_flags, int flowv_interp_,
                 int cylsymm_, int field_aligned_values_ , int _coulomb_collisions)
      : particlesPointer(_particlesPointer),
        dt(_dt),
        nR_flowV(_nR_flowV),
        nY_flowV(_nY_flowV),
        nZ_flowV(_nZ_flowV),
        flowVGridr(_flowVGridr),
        flowVGridy(_flowVGridy),
        flowVGridz(_flowVGridz),
        flowVr(_flowVr),
        flowVz(_flowVz),
        flowVt(_flowVt),
        nR_Dens(_nR_Dens),
        nZ_Dens(_nZ_Dens),
        DensGridr(_DensGridr),
        DensGridz(_DensGridz),
        ni(_ni),
        nR_Temp(_nR_Temp),
        nZ_Temp(_nZ_Temp),
        TempGridr(_TempGridr),
        TempGridz(_TempGridz),
        ti(_ti),
        te(_te),
        background_Z(_background_Z),
        background_amu(_background_amu),
        nR_Bfield(_nR_Bfield),
        nZ_Bfield(_nZ_Bfield),
        BfieldGridR(_BfieldGridR),
        BfieldGridZ(_BfieldGridZ),
        BfieldR(_BfieldR),
        BfieldZ(_BfieldZ),
        BfieldT(_BfieldT),
	gitr_flags(_gitr_flags),
        dv{0.0, 0.0, 0.0},
        state(_state),
        flowv_interp( flowv_interp_ ),
        cylsymm( cylsymm_ ),
        field_aligned_values( field_aligned_values_ ),
        coulomb_collisions(_coulomb_collisions)

        { }

CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) { 

  //bool use_ion_neutral = false;
  //if (use_ion_neutral)
  //{
  //  if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0)
  //  {
  //  gitr_precision vx = particlesPointer->vx[indx];
  //  gitr_precision vy = particlesPointer->vy[indx];
  //  gitr_precision vz = particlesPointer->vz[indx];
  //  gitr_precision charge = particlesPointer->charge[indx];
  //  gitr_precision amu = particlesPointer->amu[indx];
  //  gitr_precision mu = amu*background_amu/(amu+background_amu);
  //    double sigma = 1.5e-11;
  //    double neutral_density = 1.0e19;
  //    double v_norm = std::sqrt(vx*vx + vy*vy + vz*vz);
  //    double nu = neutral_density*sigma*v_norm;
  //    double P = 1.0 - std::exp(-dt*nu);
  //    double r1 = curand_uniform(&state[indx]);
  //    if ( r1 < P )
  //      {
  //        double psi = 2*M_PI*curand_uniform(&state[indx]);
  //        double chi = std::acos(1 - 2*curand_uniform(&state[indx]));
  //        double ux = vx;
  //        double uy = vy;
  //        double uz = vz;
  //        double u = std::sqrt(ux*ux + uy*uy + uz*uz);
  //        double u_perp = std::sqrt(ux*ux + uy*uy);
 
  //        if (u_perp == 0.0) u_perp = 1.0;
  //        if (ux == 0.0) ux = 1.0;
  //        if (uy == 0.0) uy = 1.0;
  //        if (uz == 0.0) uz = 1.0;

  //        double d_ux = ux/u_perp*uz*std::sin(chi)*std::cos(psi) - uy/u_perp*u*std::sin    (chi)*std::sin(psi) - ux*(1-std::cos(chi));
  //        double d_uy = uy/u_perp*uz*std::sin(chi)*std::cos(psi) + ux/u_perp*u*std::sin    (chi)*std::sin(psi) - uy*(1.0-std::cos(chi));
  //        double d_uz = -u_perp*std::sin(chi)*std::cos(psi) - uz*(1.0-std::cos(chi));
 
  //           particlesPointer->vx[indx] = particlesPointer->vx[indx] + mu/amu*d_ux;
  //           particlesPointer->vy[indx] = particlesPointer->vy[indx] + mu/amu*d_uy;
  //           particlesPointer->vz[indx] = particlesPointer->vz[indx] + mu/amu*d_uz;
  //      }
  //  }
  //}
  // Hard-coded option to use binary collision operator or not
  //bool use_bca = false;
  if(particlesPointer->hitWall[indx] == 0.0 && particlesPointer->charge[indx] != 0.0)
  {
    if(gitr_flags->USE_ADAPTIVE_DT)
    {
      dt = particlesPointer->dt[indx];	   
    }

  if ( coulomb_collisions == 2)
  {
    // Physical constants (should be replaced with system/better precision values)
    double ME = 9.10938356e-31;
    double MI = 1.6737236e-27;
    double Q = 1.60217662e-19;
    double EPS0 = 8.854187e-12;
    double PI = 3.141592653589793;

    // Get particle attributes
    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    gitr_precision vx = particlesPointer->vx[indx];
    gitr_precision vy = particlesPointer->vy[indx];
    gitr_precision vz = particlesPointer->vz[indx];
    gitr_precision charge = particlesPointer->charge[indx];
    gitr_precision amu = particlesPointer->amu[indx];
    gitr_precision mu = amu*background_amu/(amu+background_amu);
    gitr_precision flowVelocity[3]= {0.0};
    
    // Interpolate ion temperature
    gitr_precision ti_eV = interp2dCombined(x, y, z, nR_Temp, nZ_Temp, TempGridr, TempGridz, ti,cylsymm);
    
    // Interpolate ion density
    gitr_precision density = interp2dCombined(x, y, z, nR_Dens, nZ_Dens, DensGridr, DensGridz, ni,cylsymm);
  
    if( gitr_flags->USE_SHEATH_DENSITY )
    {
     density = density*particlesPointer->f_psi[indx];
    }

    // Calculate standard deviation for Gaussian - representing Maxwellian velocity distribution
    gitr_precision standard_deviation = 1.0/std::sqrt(2*background_amu*1.6737236e-27/(2*ti_eV*1.60217662e-19));

#ifdef __CUDACC__
    gitr_precision r1 = curand_normal(&state[indx]);
    gitr_precision r2 = curand_normal(&state[indx]);
    gitr_precision r3 = curand_normal(&state[indx]);
    gitr_precision r4 = curand_normal(&state[indx]);
    gitr_precision r5 = curand_uniform(&state[indx]);
#else
    std::normal_distribution<gitr_precision> distribution(0.0,1.0);
    gitr_precision r1 = distribution(state[indx]);
    gitr_precision r2 = distribution(state[indx]);
    gitr_precision r3 = distribution(state[indx]);
    gitr_precision r4 = distribution(state[indx]);
    std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
    gitr_precision r5 = dist(state[indx]);
#endif
    // Draw random normal background ion velocities to use in collision
    double ux_gas = standard_deviation*r1;
    double uy_gas = standard_deviation*r2;
    double uz_gas = standard_deviation*r3;

    // Interpolate flow velocity - for this problem, it is set to zero
    interp2dVector(flowVelocity,particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                        nR_flowV,nZ_flowV,
                        flowVGridr,flowVGridz,flowVr,flowVz,flowVt, cylsymm );
             
    // Shift background ion velocities by bulk flow velocities
    ux_gas = ux_gas + flowVelocity[0];
    uy_gas = uy_gas + flowVelocity[1];
    uz_gas = uz_gas + flowVelocity[2];

    // Get relative velocities and norms
    double ux = vx - ux_gas;
    double uy = vy - uy_gas;
    double uz = vz - uz_gas;
    double v_norm = std::sqrt(vx*vx + vy*vy + vz*vz);
       
    double u = std::sqrt(ux*ux + uy*uy + uz*uz);
    double u_perp = std::sqrt(ux*ux + uy*uy);

    // Calculation of the coulomb logarithm
    double lam_d = std::sqrt(EPS0*ti_eV/(density*Q));//%only one q in order to convert to J
    double lam = 12*PI*density*std::pow(lam_d,3);
    
    // Collision frequency
    double nu_0 = (1.0/std::pow(u,3.0))*std::pow(Q,4.0)*charge*charge*background_Z*background_Z*std::log(lam)*density/((mu*mu*MI*MI)*8.0*PI*EPS0*EPS0);

    // Calculate the scattering angles
    double chi_squared = dt*nu_0;
    //gitr_precision r1 = curand_uniform(&state[indx]);
    //double chi = std::sqrt(-2.0*chi_squared*std::log(r1));
    //sampling normal dist. approach
    double delta = std::sqrt(chi_squared)*r4;
    double chi = 2.0*std::atan(delta);
    double psi = 2.0*PI*r5;

    double d_ux = 0.0;
    double d_uy = 0.0;
    double d_uz = 0.0;
    // Special case where u_perp = 0
    if (u_perp == 0.0)
    {
             d_ux = u*std::sin(chi)*std::cos(psi);
             d_uy = u*std::sin(chi)*std::sin(psi);
             d_uz = -u*(1.0-std::cos(chi));

    }
    //if (u_perp == 0.0) u_perp = 1.0;
    //if (ux == 0.0) ux = 1.0;
    //if (uy == 0.0) uy = 1.0;
    //if (uz == 0.0) uz = 1.0;
    else {    
             d_ux = ux/u_perp*uz*std::sin(chi)*std::cos(psi) - uy/u_perp*u*std::sin(chi)*std::sin(psi) - ux*(1-std::cos(chi));
             d_uy = uy/u_perp*uz*std::sin(chi)*std::cos(psi) + ux/u_perp*u*std::sin(chi)*std::sin(psi) - uy*(1.0-std::cos(chi));
             d_uz = -u_perp*std::sin(chi)*std::cos(psi) - uz*(1.0-std::cos(chi));
    }
    
    if(gitr_flags->USE_ADAPTIVE_DT)
      {
        if (particlesPointer->advance[indx])
        {

             particlesPointer->vx[indx] = particlesPointer->vx[indx] + mu/amu*d_ux;
             particlesPointer->vy[indx] = particlesPointer->vy[indx] + mu/amu*d_uy;
             particlesPointer->vz[indx] = particlesPointer->vz[indx] + mu/amu*d_uz;
             particlesPointer->test[indx] = particlesPointer->vz[indx] - vz;
        }
      }
    else
    {
             particlesPointer->vx[indx] = particlesPointer->vx[indx] + mu/amu*d_ux;
             particlesPointer->vy[indx] = particlesPointer->vy[indx] + mu/amu*d_uy;
             particlesPointer->vz[indx] = particlesPointer->vz[indx] + mu/amu*d_uz;
             particlesPointer->test[indx] = particlesPointer->vz[indx] - vz;
    }
    
  }
  else if (coulomb_collisions == 1)
  {
     
    gitr_precision pi = 3.14159265;   
    //gitr_precision k_boltz = 1.38e-23*11604/1.66e-27;
    gitr_precision T_background = 0.0;
    gitr_precision nu_friction = 0.0;
    gitr_precision nu_deflection = 0.0;
    gitr_precision nu_parallel = 0.0;
    gitr_precision nu_energy = 0.0;
    gitr_precision flowVelocity[3]= {0.0};
    gitr_precision relativeVelocity[3] = {0.0};
    gitr_precision velocityCollisions[3]= {0.0};	
    gitr_precision velocityRelativeNorm;	
    gitr_precision parallel_direction[3] = {0.0};
    gitr_precision perp_direction1[3] = {0.0};
    gitr_precision perp_direction2[3] = {0.0};
    
    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    gitr_precision vx = particlesPointer->vx[indx];
    gitr_precision vy = particlesPointer->vy[indx];
    gitr_precision vz = particlesPointer->vz[indx];

    if( flowv_interp == 3 )
    {
    interp3dVector (&flowVelocity[0], particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],nR_flowV,nY_flowV,nZ_flowV,
                flowVGridr,flowVGridy,flowVGridz,flowVr,flowVz,flowVt);
    }
    else if( flowv_interp < 3 )
    {
    if( field_aligned_values > 0 )
    {
    interpFieldAlignedVector(&flowVelocity[0],
                                 particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                                 nR_flowV,nZ_flowV,
                                 flowVGridr,flowVGridz,flowVr,
                                 flowVz,flowVt,nR_Bfield,nZ_Bfield,
                                 BfieldGridR,BfieldGridZ,BfieldR,
                                 BfieldZ,BfieldT, cylsymm );
    }
    else
    {
    interp2dVector(flowVelocity,particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                        nR_flowV,nZ_flowV,
                        flowVGridr,flowVGridz,flowVr,flowVz,flowVt, cylsymm );
    }
    }

    relativeVelocity[0] = vx - flowVelocity[0];
    relativeVelocity[1] = vy - flowVelocity[1];
    relativeVelocity[2] = vz - flowVelocity[2];
    velocityRelativeNorm = vectorNorm(relativeVelocity);

#ifdef __CUDACC__
    gitr_precision n1 = curand_normal(&state[indx]);
    gitr_precision n2 = curand_normal(&state[indx]);
    gitr_precision r1 = curand_uniform(&state[indx]);
    gitr_precision r2 = curand_uniform(&state[indx]);
    gitr_precision r3 = curand_uniform(&state[indx]);
    gitr_precision xsi = curand_uniform(&state[indx]);
#else
    std::normal_distribution<gitr_precision> distribution(0.0,1.0);
    std::uniform_real_distribution<gitr_precision> dist(0.0, 1.0);
    gitr_precision n1 = distribution(state[indx]);
    gitr_precision n2 = distribution(state[indx]);
    gitr_precision r1 = dist(state[indx]);
    gitr_precision r2 = dist(state[indx]);
    gitr_precision r3 = dist(state[indx]);
    gitr_precision xsi = dist(state[indx]);
#endif

    getSlowDownFrequencies(nu_friction, nu_deflection, nu_parallel, nu_energy,
                             x, y, z,
                             vx, vy, vz,
                             particlesPointer->charge[indx], particlesPointer->amu[indx],
                             nR_flowV, nZ_flowV, flowVGridr,
                             flowVGridz, flowVr,
                             flowVz, flowVt,
                             nR_Dens, nZ_Dens, DensGridr,
                             DensGridz, ni, nR_Temp, nZ_Temp,
                             TempGridr, TempGridz, ti, te, background_Z, background_amu,
                             nR_Bfield,
                             nZ_Bfield,
                             BfieldGridR,
                             BfieldGridZ,
                             BfieldR,
                             BfieldZ,
                             BfieldT, T_background, flowv_interp, cylsymm,
                             field_aligned_values,gitr_flags->USE_SHEATH_DENSITY,particlesPointer->f_psi[indx]  );

    getSlowDownDirections2(parallel_direction, perp_direction1, perp_direction2,
                            relativeVelocity[0] , relativeVelocity[1] , relativeVelocity[2] );
      
    gitr_precision ti_eV = interp2dCombined( x, y, z, nR_Temp, nZ_Temp, TempGridr, TempGridz, ti,
                                             cylsymm );

    gitr_precision density = interp2dCombined( x, y, z, nR_Dens, nZ_Dens, DensGridr, 
                                               DensGridz, ni, cylsymm );
    
    //printf("speed %f temp %f density %4.3e charge %2.2f nu_slowing_down %4.3e \n",velocityRelativeNorm,ti_eV,density,particlesPointer->charge[indx], nu_friction);
    //printf("nu_deflection %4.3e nu_parallel %4.3e nu_energy %4.3e \n",nu_deflection, nu_parallel, nu_energy);
    if(nu_parallel <=0.0) nu_parallel = 0.0;
    gitr_precision coeff_par = n1 * std::sqrt(2.0*nu_parallel * dt);
    gitr_precision cosXsi = cos(2.0 * pi * xsi) - 0.0028;
    if(cosXsi > 1.0) cosXsi = 1.0;
    gitr_precision sinXsi = sin(2.0 * pi * xsi);
    if(nu_deflection <=0.0) nu_deflection = 0.0;
    gitr_precision coeff_perp1 = cosXsi * std::sqrt(nu_deflection * dt*0.5);
    gitr_precision coeff_perp2 = sinXsi * std::sqrt(nu_deflection * dt*0.5);
      
    gitr_precision nuEdt = nu_energy * dt;
    if (nuEdt < -1.0) nuEdt = -1.0;
      
    gitr_precision vx_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[0] + std::abs(n2)*(coeff_perp1 * perp_direction1[0] + coeff_perp2 * perp_direction2[0])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[0];
    gitr_precision vy_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[1] + std::abs(n2)*(coeff_perp1 * perp_direction1[1] + coeff_perp2 * perp_direction2[1])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[1];
    gitr_precision vz_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[2] + std::abs(n2)*(coeff_perp1 * perp_direction1[2] + coeff_perp2 * perp_direction2[2])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[2];
      
    if(gitr_flags->USE_ADAPTIVE_DT)
      {
        if (particlesPointer->advance[indx])
        {
          particlesPointer->vx[indx] = vx_relative + flowVelocity[0]; 
          particlesPointer->vy[indx] = vy_relative + flowVelocity[1]; 
          particlesPointer->vz[indx] = vz_relative + flowVelocity[2];
        }
      }
      else
      {
        particlesPointer->vx[indx] = vx_relative + flowVelocity[0]; 
        particlesPointer->vy[indx] = vy_relative + flowVelocity[1]; 
        particlesPointer->vz[indx] = vz_relative + flowVelocity[2];
      }

      this->dv[0] = velocityCollisions[0];
      this->dv[1] = velocityCollisions[1];
      this->dv[2] = velocityCollisions[2];
    }
    else if (coulomb_collisions == 3)
    {
    // Physical constants (should be replaced with system/better precision values)
    double ME = 9.10938356e-31;
    double MI = 1.6737236e-27;
    double Q = 1.60217662e-19;
    double EPS0 = 8.854187e-12;
    double PI = 3.141592653589793;

    // Get particle attributes
    gitr_precision x = particlesPointer->xprevious[indx];
    gitr_precision y = particlesPointer->yprevious[indx];
    gitr_precision z = particlesPointer->zprevious[indx];
    gitr_precision vx = particlesPointer->vx[indx];
    gitr_precision vy = particlesPointer->vy[indx];
    gitr_precision vz = particlesPointer->vz[indx];
    gitr_precision charge = particlesPointer->charge[indx];
    gitr_precision amu = particlesPointer->amu[indx];
    gitr_precision mu = amu*background_amu/(amu+background_amu);
    gitr_precision v_flow[3]= {0.0};
    gitr_precision b[3]= {0.0};
    
    gitr_precision nu_s = 0.0;
    gitr_precision nu_d = 0.0;
    gitr_precision nu_par = 0.0;
    gitr_precision nu_E = 0.0;
    
    gitr_precision T_background = 0.0;
   // // Interpolate ion temperature
   // gitr_precision ti_eV = interp2dCombined(x, y, z, nR_Temp, nZ_Temp, TempGridr, TempGridz, ti,cylsymm);
   // 
   // // Interpolate ion density
   // gitr_precision density = interp2dCombined(x, y, z, nR_Dens, nZ_Dens, DensGridr, DensGridz, ni,cylsymm);
  
   // if( use_sheath_density )
   // {
   //  density = density*particlesPointer->f_psi[indx];
   // }
    
    interp2dVector(v_flow,particlesPointer->xprevious[indx],particlesPointer->yprevious[indx],particlesPointer->zprevious[indx],
                        nR_flowV,nZ_flowV,
                        flowVGridr,flowVGridz,flowVr,flowVz,flowVt, cylsymm );
    
    gitr_precision wx = vx - v_flow[0];
    gitr_precision wy = vy - v_flow[1];
    gitr_precision wz = vz - v_flow[2];

    gitr_precision vnorm = std::sqrt(v_flow[0]*v_flow[0] + v_flow[1]*v_flow[1] + v_flow[2]*v_flow[2]);

    if (vnorm < 1.0e-12)
    {
      b[0] = 0;
      b[1] = 1;
      b[2] = 0;
    }
    else
    {
      b[0] = v_flow[0]/vnorm;
      b[1] = v_flow[1]/vnorm;
      b[2] = v_flow[2]/vnorm;
    }

    gitr_precision w_norm = std::sqrt(wx*wx + wy*wy + wz*wz);
    gitr_precision w1x = wx/w_norm;
    gitr_precision w1y = wy/w_norm;
    gitr_precision w1z = wz/w_norm;

    if (w1x == b[0] && w1y == b[1] && w1z == b[2])
    {
        wz = wz + 1e-2;
    }

    gitr_precision w_norm_squared = wx*wx + wy*wy + wz*wz;
    w_norm = std::sqrt(w_norm_squared);

    w1x = wx/w_norm;
    w1y = wy/w_norm;
    w1z = wz/w_norm;

    gitr_precision w1_dot_b = b[0]*w1x + b[1]*w1y +b[2]*w1z;
    
    gitr_precision sqrt_term = sqrt(1-w1_dot_b*w1_dot_b);

    gitr_precision w2x = 1.0/sqrt_term*(w1_dot_b*w1x - b[0]);
    gitr_precision w2y = 1.0/sqrt_term*(w1_dot_b*w1y - b[1]);
    gitr_precision w2z = 1.0/sqrt_term*(w1_dot_b*w1z - b[2]);

    gitr_precision w3x = 1.0/sqrt_term*(b[2]*w1y - b[1]*w1z);
    gitr_precision w3y = 1.0/sqrt_term*(b[0]*w1z - b[2]*w1x);
    gitr_precision w3z = 1.0/sqrt_term*(b[1]*w1x - b[0]*w1y);

//[nu_s nu_d nu_par nu_E] = getFrequencies(w_norm,T,m,mD,nD,z,zD);
    getSlowDownFrequencies(nu_s, nu_d, nu_par, nu_E,
                             x, y, z,
                             vx, vy, vz,
                             particlesPointer->charge[indx], particlesPointer->amu[indx],
                             nR_flowV, nZ_flowV, flowVGridr,
                             flowVGridz, flowVr,
                             flowVz, flowVt,
                             nR_Dens, nZ_Dens, DensGridr,
                             DensGridz, ni, nR_Temp, nZ_Temp,
                             TempGridr, TempGridz, ti, te, background_Z, background_amu,
                             nR_Bfield,
                             nZ_Bfield,
                             BfieldGridR,
                             BfieldGridZ,
                             BfieldR,
                             BfieldZ,
                             BfieldT, T_background, flowv_interp, cylsymm,
                             field_aligned_values,gitr_flags->USE_ADAPTIVE_DT,particlesPointer->f_psi[indx]  );


//gamma1 = normrnd(0,1,nP,1);
//gamma2 = normrnd(0,1,nP,1);
//gamma3 = normrnd(0,1,nP,1);
#ifdef __CUDACC__
    gitr_precision gamma1 = curand_normal(&state[indx]);
    gitr_precision gamma2 = curand_normal(&state[indx]);
    gitr_precision gamma3 = curand_normal(&state[indx]);
#else
    std::normal_distribution<gitr_precision> distribution(0.0,1.0);
    gitr_precision gamma1 = distribution(state[indx]);
    gitr_precision gamma2 = distribution(state[indx]);
    gitr_precision gamma3 = distribution(state[indx]);
#endif
    gitr_precision term1 = dt*(-nu_s*w_norm);
    gitr_precision term2 = std::sqrt(dt*w_norm_squared*nu_par)*gamma1;
    gitr_precision term3 = std::sqrt(0.5*dt*w_norm_squared*nu_d)*gamma2;
    gitr_precision term4 = std::sqrt(0.5*dt*w_norm_squared*nu_d)*gamma3;

gitr_precision dvx = w1x*term1 + w1x*term2 + w2x*term3 + w3x*term4;
gitr_precision dvy = w1y*term1 + w1y*term2 + w2y*term3 + w3y*term4;
gitr_precision dvz = w1z*term1 + w1z*term2 + w2z*term3 + w3z*term4;

    if(gitr_flags->USE_ADAPTIVE_DT)
      {
        if (particlesPointer->advance[indx])
        {
        particlesPointer->vx[indx] = vx + dvx; 
        particlesPointer->vy[indx] = vy + dvy; 
        particlesPointer->vz[indx] = vz + dvz;
        }
      }
      else
      {
        particlesPointer->vx[indx] = vx + dvx; 
        particlesPointer->vy[indx] = vy + dvy; 
        particlesPointer->vz[indx] = vz + dvz;
      }

    }
  }
  }
};

#endif
