//------------------------------------------------------------------------------
// GITR: thermalForce.h
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - August 2023 by Diaw
//
// Description:
//     This header file defines thermal force
//
// Models Included:
//     1. Stangeby Model:
//        - Reference: 
//            - P.C. Stangeby "The Plasma Boundary of Magnetic Fusion Devices", (see page 299-300)
//             F_e =beta grad T_e where beta is function of reduced mass and Z
//             F_i =alpha grad T_i  where alpha = 0.71 Z^2
//
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------


#ifndef _THERMAL_
#define _THERMAL_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "Particles.h"
#include <cmath>
#include "constants.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct thermalForce { 
    Flags* flags;
    Particles* p;
    const gitr_precision dt;
    gitr_precision background_amu;
    int nR_gradT, nZ_gradT, nR_Bfield, nZ_Bfield, cylsymm;
    
    gitr_precision* gradTGridr;
    gitr_precision* gradTGridz;
    gitr_precision* gradTiR;
    gitr_precision* gradTiZ;
    gitr_precision* gradTiT;
    gitr_precision* gradTeR;
    gitr_precision* gradTeZ;
    gitr_precision* gradTeT;
    gitr_precision* BfieldGridRDevicePointer;
    gitr_precision* BfieldGridZDevicePointer;
    gitr_precision* BfieldRDevicePointer;
    gitr_precision* BfieldZDevicePointer;
    gitr_precision* BfieldTDevicePointer;

 thermalForce(
        Flags* _flags,
        Particles* _p,
        gitr_precision _dt,
        gitr_precision _background_amu,
        int _nR_gradT,
        int _nZ_gradT,
        gitr_precision* _gradTGridr,
        gitr_precision* _gradTGridz,
        gitr_precision* _gradTiR,
        gitr_precision* _gradTiZ,
        gitr_precision* _gradTiT,
        gitr_precision* _gradTeR,
        gitr_precision* _gradTeZ,
        gitr_precision* _gradTeT,
        int _nR_Bfield,
        int _nZ_Bfield,
        gitr_precision* _BfieldGridRDevicePointer,
        gitr_precision* _BfieldGridZDevicePointer,
        gitr_precision* _BfieldRDevicePointer,
        gitr_precision* _BfieldZDevicePointer,
        gitr_precision* _BfieldTDevicePointer,
        int _cylsymm
    )
    : flags(_flags),
      p(_p),
      dt(_dt),
      background_amu(_background_amu),
      nR_gradT(_nR_gradT),
      nZ_gradT(_nZ_gradT),
      gradTGridr(_gradTGridr),
      gradTGridz(_gradTGridz),
      gradTiR(_gradTiR),
      gradTiZ(_gradTiZ),
      gradTiT(_gradTiT),
      gradTeR(_gradTeR),
      gradTeZ(_gradTeZ),
      gradTeT(_gradTeT),
      nR_Bfield(_nR_Bfield),
      nZ_Bfield(_nZ_Bfield),
      BfieldGridRDevicePointer(_BfieldGridRDevicePointer),
      BfieldGridZDevicePointer(_BfieldGridZDevicePointer),
      BfieldRDevicePointer(_BfieldRDevicePointer),
      BfieldZDevicePointer(_BfieldZDevicePointer),
      BfieldTDevicePointer(_BfieldTDevicePointer),
      cylsymm(_cylsymm) 
    {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx) { 
    if (p->hitWall[indx] == 0.0 && p->charge[indx] > 0.0) {
        gitr_precision alpha, beta, mu;
        gitr_precision gradTe[3] = {0.0}, gradTi[3] = {0.0}, B[3] = {0.0}, B_unit[3] = {0.0};
        gitr_precision Bmag, gradTiPar, dv_ITG[3] = {}, dv_ETG[3] = {}, dt_step = dt;

        if (flags->USE_ADAPTIVE_DT) {
            dt_step = p->advance[indx] ? p->dt[indx] : 0.0;
        }

        // Interpolation functions to get gradTe and gradTi
        interp2dVector(gradTi, p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT, gradTGridr, gradTGridz, gradTiR, gradTiZ, gradTiT, cylsymm);
        interp2dVector(gradTe, p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_gradT, nZ_gradT,gradTGridr, gradTGridz, gradTeR, gradTeZ, gradTeT, cylsymm);
        
        // Calculate alpha and beta
        mu = p->amu[indx] / (background_amu + p->amu[indx]);
        alpha = p->charge[indx] * p->charge[indx] * 0.71;
        beta = 3.0 * (mu + 5.0 * std::sqrt(2.0) * p->charge[indx] * p->charge[indx] * (1.1 * std::pow(mu, 2.5) - 0.35 * std::pow(mu, 1.5)) - 1) / (2.6 - 2 * mu + 5.4 * std::pow(mu, 2));

        interp2dVector(B, p->xprevious[indx], p->yprevious[indx], p->zprevious[indx], nR_Bfield, nZ_Bfield, BfieldGridRDevicePointer, BfieldGridZDevicePointer, BfieldRDevicePointer, BfieldZDevicePointer, BfieldTDevicePointer, cylsymm);    
        Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
        
        for(int i = 0; i < 3; i++) {
            B_unit[i] = B[i] / Bmag;
            dv_ETG[i] = gitr_constants::e * dt_step / (p->amu[indx] * gitr_constants::m_p) * (alpha * gradTe[i]) * B_unit[i];
            dv_ITG[i] = gitr_constants::e * dt_step / (p->amu[indx] * gitr_constants::m_p) * (beta * gradTi[i]) * B_unit[i];
        }
   
        p->vx[indx] += dv_ITG[0];
        p->vy[indx] += dv_ITG[1];
        p->vz[indx] += dv_ITG[2];
    }
}
};
#endif
