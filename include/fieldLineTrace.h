#ifndef _FIELDTRACE_
#define _FIELDTRACE_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif
#include "Particles.h"
#include "Boundary.h"
#include "boris.h"
#include <cmath>

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct field_line_trace {
    gitr_precision BfieldFactor; 
    Particles *particles;
    gitr_precision dr;
    Boundary *boundaries;
    int nLines;
    int nR_Lc;
    int nZ_Lc;
    gitr_precision* gridRLc;
    gitr_precision* gridZLc;
    gitr_precision* Lc;
            int nR_Bfield;
            int nZ_Bfield;
            gitr_precision * BfieldGridR;
            gitr_precision * BfieldGridZ;
            gitr_precision * BfieldR;
            gitr_precision * BfieldZ;
            gitr_precision * BfieldT;
    int cylsymm;
            
    field_line_trace(gitr_precision _BfieldFactor,Particles* _particles,gitr_precision _dr,Boundary* _boundaries,int _nLines, int _nR_Lc, int _nZ_Lc, 
            gitr_precision* _gridRLc, gitr_precision* _gridZLc, gitr_precision* _Lc,
            int _nR_Bfield, int _nZ_Bfield,
            gitr_precision * _BfieldGridR,
            gitr_precision * _BfieldGridZ,
            gitr_precision * _BfieldR,
            gitr_precision * _BfieldZ,
            gitr_precision * _BfieldT,
            int cylsymm_ )
        
            : BfieldFactor(_BfieldFactor),particles(_particles),dr(_dr),boundaries(_boundaries),nLines(_nLines),
        nR_Lc(_nR_Lc),nZ_Lc(_nZ_Lc),
        gridRLc(_gridRLc), gridZLc(_gridZLc),Lc(_Lc),
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridR(_BfieldGridR), BfieldGridZ(_BfieldGridZ),
    BfieldR(_BfieldR), BfieldZ(_BfieldZ), BfieldT(_BfieldT), cylsymm( cylsymm_ ) {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx) const { 
    gitr_precision B[3] = {0.0f,0.0f,0.0f};
    gitr_precision Bnorm[3] = {0.0f,0.0f,0.0f};
    gitr_precision Bmag = 0.0f;
    gitr_precision particleDistance = 0.0f;
    gitr_precision k1[3] = {0.0,0.0,0.0};
    gitr_precision k2[3] = {0.0,0.0,0.0};
    gitr_precision k3[3] = {0.0,0.0,0.0};
    gitr_precision k4[3] = {0.0,0.0,0.0};
    gitr_precision x0 = 0.0f;
    gitr_precision y0 = 0.0f;
    gitr_precision z0 = 0.0f;
    gitr_precision x = 0.0f;
    gitr_precision y = 0.0f;
    gitr_precision z = 0.0f;

    gitr_precision dr_fac = BfieldFactor*dr;

if(particles->hitWall[indx] == 0.0)
{
    x0 = particles->x[indx];
    y0 = particles->y[indx];
    z0 = particles->z[indx];

    interp2dVector(&B[0],x0, y0,z0,
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT, cylsymm );
    //std::cout << "Bfield interp " << B[0] << " " << B[1] << " " << B[2] << std::endl;
    vectorNormalize(B,B);
    //Bmag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    //Bnorm[0] = B[0]/Bmag;
    //Bnorm[1] = B[1]/Bmag;
    //Bnorm[2] = B[2]/Bmag;
    //std::cout << "Bfield interp " << B[0] << " " << B[1] << " " << B[2] << std::endl;
    vectorScalarMult(dr_fac,B,k1);

    interp2dVector(&B[0],x0+0.5*k1[0],y0+0.5*k1[1],z0+0.5*k1[2],
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT, cylsymm );
    vectorNormalize(B,B);

    vectorScalarMult(dr_fac,B,k2);
    interp2dVector(&B[0],x0+0.5*k2[0],y0+0.5*k2[1],z0+0.5*k2[2],
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT, cylsymm );
    vectorNormalize(B,B);

    vectorScalarMult(dr_fac,B,k3);

    interp2dVector(&B[0],x0+k3[0],y0+k3[1],z0+k3[2],
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT, cylsymm );
    vectorNormalize(B,B);

    vectorScalarMult(dr_fac,B,k4);
    x = x0+k1[0]/6.0+k2[0]/3.0+k3[0]/3.0+k4[0]/6.0;
    y = y0+k1[1]/6.0+k2[1]/3.0+k3[1]/3.0+k4[1]/6.0;
    z = z0+k1[2]/6.0+k2[2]/3.0+k3[2]/3.0+k4[2]/6.0;
    particles->x[indx] = x; 
    particles->y[indx] = y;
    particles->z[indx] = z;
    //particles->x[indx] = particles->xprevious[indx] + BfieldFactor*dr*Bnorm[0];
    //particles->y[indx] = particles->yprevious[indx] + BfieldFactor*dr*Bnorm[1];
    //particles->z[indx] = particles->zprevious[indx] + BfieldFactor*dr*Bnorm[2];
    particles->distanceTraveled[indx] = particles->distanceTraveled[indx] + dr; 
}
else if(particles->hitWall[indx] == 1.0)
{
    particleDistance = std::sqrt((particles->x[indx] - particles->xprevious[indx])*
                             (particles->x[indx] - particles->xprevious[indx]) 
                           + (particles->y[indx] - particles->yprevious[indx])*
                             (particles->y[indx] - particles->yprevious[indx])
                           + (particles->z[indx] - particles->zprevious[indx])*
                             (particles->z[indx] - particles->zprevious[indx]));  
            
    particles->distanceTraveled[indx] = particles->distanceTraveled[indx]
                                        + particleDistance; 
    particles->hitWall[indx] = particles->hitWall[indx]+1.0;        
}
  
       }
     
};

#endif
