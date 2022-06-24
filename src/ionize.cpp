#include "ionize.h"
#include "interpRateCoeff.hpp"

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
          int _nTemperaturesIonize,
          int _nDensitiesIonize,
          gitr_precision *_gridTemperature_Ionization,
          gitr_precision *_gridDensity_Ionization,
          gitr_precision *_rateCoeff_Ionization,
          gitr_precision *_random_uniform_number,
          int cylsymm_ )
      :

        flags(_flags), particlesPointer(_particlesPointer), nR_Dens(_nR_Dens),
        nZ_Dens(_nZ_Dens), DensGridr(_DensGridr), DensGridz(_DensGridz),
        ne(_ne), nR_Temp(_nR_Temp), nZ_Temp(_nZ_Temp), TempGridr(_TempGridr),
        TempGridz(_TempGridz), te(_te),
        nTemperaturesIonize(_nTemperaturesIonize),
        nDensitiesIonize(_nDensitiesIonize),
        gridDensity_Ionization(_gridDensity_Ionization),
        gridTemperature_Ionization(_gridTemperature_Ionization),
        rateCoeff_Ionization(_rateCoeff_Ionization),
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
    if (flags->USE_ADAPTIVE_DT)
    {
      dt = particlesPointer->dt[indx];
    }
    
    tion = interpRateCoeff2d(
          particlesPointer->charge[indx], particlesPointer->x[indx],
          particlesPointer->y[indx], particlesPointer->z[indx], nR_Temp,
          nZ_Temp, TempGridr, TempGridz, te, DensGridr, DensGridz, ne,
          nTemperaturesIonize, nDensitiesIonize, gridTemperature_Ionization,
          gridDensity_Ionization, rateCoeff_Ionization, cylsymm );
    
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

/* explicit template instantiations */
#if USE_CUDA
template struct ionize< curandState  >;
#else
template struct ionize< std::mt19937 >;
#endif
