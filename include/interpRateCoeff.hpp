#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "interp2d.hpp"
#include <thrust/device_vector.h>
#include <vector>
#include <cmath>

using namespace std;

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

CUDA_CALLABLE_MEMBER
gitr_precision rateCoeffInterp(int charge, gitr_precision te, gitr_precision ne,int nT, int nD, gitr_precision* rateGrid_Tempp,gitr_precision* rateGrid_Densp,gitr_precision* Ratesp);

CUDA_CALLABLE_MEMBER
gitr_precision interpRateCoeff2d ( int charge, gitr_precision x, gitr_precision y, gitr_precision z,int nx, int nz, gitr_precision* tempGridxp,
       gitr_precision* tempGridzp, gitr_precision* Tempp,
       gitr_precision* densGridxp,gitr_precision* densGridzp,gitr_precision* Densp,int nT_Rates, int nD_Rates,
       gitr_precision* rateGrid_Temp,gitr_precision* rateGrid_Dens,gitr_precision* Rates,
       int cylsymm );
