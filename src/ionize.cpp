//------------------------------------------------------------------------------
// GITR: ionize.cpp
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - August 2023 by Diaw
//
// Description:
//     Compute rates and performs ioniziation. 
//
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------



#include "ionize.h"
#include "interpRateCoeff.hpp"
#include "processIonizationRecombination.h"
#include "interpolator2.h"


#if USE_CUDA

/* I think these "get_rand" functions are never actually called */
CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(curandState *state,int indx)
{
  curandState localState = state[indx];
  float r1 = curand_uniform(&localState);
  state[indx] = localState;
  return r1;
}

#else

CUDA_CALLABLE_MEMBER_HOST CUDA_CALLABLE_MEMBER_DEVICE
float get_rand(std::mt19937 *state,int indx)
{
  std::uniform_real_distribution<float> dist(0.0, 1.0);
  std::mt19937 this_state = state[indx];
  float r1 = dist(this_state);
  state[indx] = this_state;
  return r1;
}

#endif

#if __CUDA_ARCH__
CUDA_CALLABLE_MEMBER_DEVICE
double get_rand_double(curandState *state,int indx)
{
  curandState localState = state[indx];
  double r1 = curand_uniform_double(&localState);
  state[indx] = localState;
  return r1;
}
#endif

double get_rand_double(std::mt19937 *state,int indx)
{
  std::uniform_real_distribution<double> dist(-1.0, 0.0);
  std::mt19937 this_state = state[indx];
  double r1 = -dist(this_state);
  state[indx] = this_state;
  return r1;
}

inline std::vector<std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>>> ionizationDataCache;


template< typename T >
ionize< T >::ionize(Flags *_flags,
          Particles *_particlesPointer,
          gitr_precision _dt,
          T *_state,
          int _nR_Dens,
          int _nZ_Dens,
          gitr_precision *_DensGridr,
          gitr_precision *_DensGridz,
          gitr_precision *_ne,
          int _nR_Temp,
          int _nZ_Temp,
          gitr_precision *_TempGridr,
          gitr_precision *_TempGridz,
          gitr_precision *_te,
          gitr_precision *_random_uniform_number,
          int cylsymm_ )
      :

        flags(_flags), particlesPointer(_particlesPointer), nR_Dens(_nR_Dens),
        nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),
        ne(_ne), nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr),
        TempGridz(_TempGridz), te(_te),
        dt(_dt),
        state(_state),random_uniform_number{_random_uniform_number},
        cylsymm( cylsymm_ )
{ }


template< typename T >
CUDA_CALLABLE_MEMBER_DEVICE
void ionize< T >::operator()(std::size_t indx)
{
  if (flags->USE_IONIZATION)
  {
    // Get ionization data from ionizationDataCache or process it
    std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> ionizationResult;
    if(particlesPointer->Z[indx] < ionizationDataCache.size() && ionizationDataCache[particlesPointer->Z[indx]] != std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>>()) {
        ionizationResult = ionizationDataCache[particlesPointer->Z[indx]];
    } else {
        ionizationResult = process_rates(IONIZATION, particlesPointer->Z[indx]);
        // Resize the cache vector if needed
        if(particlesPointer->Z[indx] >= ionizationDataCache.size()) {
            ionizationDataCache.resize(particlesPointer->Z[indx] + 1);
        }
        ionizationDataCache[particlesPointer->Z[indx]] = ionizationResult;
    }
    
    size_t nTemperaturesIonize = std::get<0>(ionizationResult);
    size_t nDensitiesIonize = std::get<1>(ionizationResult);
    size_t nChargeStatesIonize = std::get<2>(ionizationResult);
    const std::vector<double>& gridTemperature_Ionization = std::get<3>(ionizationResult);
    const std::vector<double>& gridDensity_Ionization = std::get<4>(ionizationResult);
    const std::vector<double>& rateCoeff_Ionization = std::get<5>(ionizationResult);


    if (flags->USE_ADAPTIVE_DT)
    {
      dt = particlesPointer->dt[indx];
    }

    // Get local temperature and density
    gitr_precision tlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx], nR_Temp, nZ_Temp, TempGridr, TempGridz, te);
    gitr_precision nlocal = interp2dCombined(particlesPointer->x[indx], particlesPointer->y[indx], particlesPointer->z[indx], nR_Dens, nZ_Dens, DensGridr, DensGridz, ne);
    gitr_precision RClocal = rateCoeffInterp(particlesPointer->charge[indx], tlocal, nlocal, nTemperaturesIonize, nDensitiesIonize, gridTemperature_Ionization, gridDensity_Ionization, rateCoeff_Ionization);

    tion = 1.0 / (RClocal * nlocal);
    gitr_precision P = exp(-dt / tion);
    gitr_precision P1 = 1.0 - P;
    gitr_precision r1 = get_rand_double(state,indx);
    random_uniform_number[0] = r1;
    
    if (flags->USE_ADAPTIVE_DT) 
    {
      if (particlesPointer->hitWall[indx] == 0.0 && particlesPointer->advance[indx])
      {
        if (r1 <= P1) 
        {
          particlesPointer->charge[indx] = particlesPointer->charge[indx] + 1;
        }
      }
    }
    else
    {
      if (particlesPointer->hitWall[indx] == 0.0)
      {
        if (r1 <= P1)
        {
          particlesPointer->charge[indx] = particlesPointer->charge[indx] + 1;
        }
      }
    }
  }
}


/* This function interpolates the rate coefficient from the rate coefficient grid */
gitr_precision rateCoeffInterp(int charge, gitr_precision te, gitr_precision ne, int nT, int nD, const std::vector<double>& rateGrid_Tempp, const std::vector<double>& rateGrid_Densp, const std::vector<double>& Ratesp)
{
    int indT = 0;
    int indN = 0;
    gitr_precision logT = std::log10(te);
    gitr_precision logn = std::log10(ne);
    gitr_precision d_T = rateGrid_Tempp[1] - rateGrid_Tempp[0];
    gitr_precision d_n = rateGrid_Densp[1] - rateGrid_Densp[0];
    indT = std::floor((logT - rateGrid_Tempp[0]) / d_T);
    indN = std::floor((logn - rateGrid_Densp[0]) / d_n);
    if (indT < 0 || indT > nT - 2) {
        indT = 0;
    }
    if (indN < 0 || indN > nD - 2) {
        indN = 0;
    }
    if (charge > 74 - 1) {
        charge = 0;
    }
    gitr_precision aT = std::pow(10.0, rateGrid_Tempp[indT + 1]) - te;
    gitr_precision bT = te - std::pow(10.0, rateGrid_Tempp[indT]);
    gitr_precision abT = aT + bT;
    gitr_precision aN = std::pow(10.0, rateGrid_Densp[indN + 1]) - ne;
    gitr_precision bN = ne - std::pow(10.0, rateGrid_Densp[indN]);
    gitr_precision abN = aN + bN;
    gitr_precision fx_z1 = (aN * std::pow(10.0, Ratesp[charge * nT * nD + indT * nD + indN]) +
        bN * std::pow(10.0, Ratesp[charge * nT * nD + indT * nD + indN + 1])) / abN;
    gitr_precision fx_z2 = (aN * std::pow(10.0, Ratesp[charge * nT * nD + (indT + 1) * nD + indN]) +
        bN * std::pow(10.0, Ratesp[charge * nT * nD + (indT + 1) * nD + indN + 1])) / abN;
    gitr_precision fxz = (aT * fx_z1 + bT * fx_z2) / abT;
    return fxz;
}



gitr_precision interp2dCombined ( gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz,
    gitr_precision* gridx,gitr_precision* gridz,gitr_precision* data ) {
    
    gitr_precision fxz = 0.0;
    gitr_precision fx_z1 = 0.0;
    gitr_precision fx_z2 = 0.0; 
    if(nx*nz == 1)
    {
        fxz = data[0];
    }
    else{
    gitr_precision dim1;
    dim1 = x;
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
    }
    }

    return fxz;
}


/* explicit template instantiations */
#if USE_CUDA
template struct ionize< curandState  >;
#else
template struct ionize< std::mt19937 >;
#endif
