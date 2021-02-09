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

struct field_line_trace {
    float BfieldFactor; 
    Particles *particles;
    float dr;
    Boundary *boundaries;
    int nLines;
    int nR_Lc;
    int nZ_Lc;
    float* gridRLc;
    float* gridZLc;
    float* Lc;
            int nR_Bfield;
            int nZ_Bfield;
            float * BfieldGridR;
            float * BfieldGridZ;
            float * BfieldR;
            float * BfieldZ;
            float * BfieldT;
            
    field_line_trace(float _BfieldFactor,Particles* _particles,float _dr,Boundary* _boundaries,int _nLines, int _nR_Lc, int _nZ_Lc, 
            float* _gridRLc, float* _gridZLc, float* _Lc,
            int _nR_Bfield, int _nZ_Bfield,
            float * _BfieldGridR,
            float * _BfieldGridZ,
            float * _BfieldR,
            float * _BfieldZ,
            float * _BfieldT)
        
            : BfieldFactor(_BfieldFactor),particles(_particles),dr(_dr),boundaries(_boundaries),nLines(_nLines),
        nR_Lc(_nR_Lc),nZ_Lc(_nZ_Lc),
        gridRLc(_gridRLc), gridZLc(_gridZLc),Lc(_Lc),
             nR_Bfield(_nR_Bfield), nZ_Bfield(_nZ_Bfield), BfieldGridR(_BfieldGridR), BfieldGridZ(_BfieldGridZ),
    BfieldR(_BfieldR), BfieldZ(_BfieldZ), BfieldT(_BfieldT) {}

CUDA_CALLABLE_MEMBER    
void operator()(std::size_t indx) const { 
    float B[3] = {0.0f,0.0f,0.0f};
    float Bnorm[3] = {0.0f,0.0f,0.0f};
    float Bmag = 0.0f;
    float particleDistance = 0.0f;
    float k1[3] = {0.0,0.0,0.0};
    float k2[3] = {0.0,0.0,0.0};
    float k3[3] = {0.0,0.0,0.0};
    float k4[3] = {0.0,0.0,0.0};
    float x0 = 0.0f;
    float y0 = 0.0f;
    float z0 = 0.0f;
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    float dr_fac = BfieldFactor*dr;

if(particles->hitWall[indx] == 0.0)
{
    x0 = particles->x[indx];
    y0 = particles->y[indx];
    z0 = particles->z[indx];

    interp2dVector(&B[0],x0, y0,z0,
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT);
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
            BfieldR,BfieldZ,BfieldT);
    vectorNormalize(B,B);

    vectorScalarMult(dr_fac,B,k2);
    interp2dVector(&B[0],x0+0.5*k2[0],y0+0.5*k2[1],z0+0.5*k2[2],
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT);
    vectorNormalize(B,B);

    vectorScalarMult(dr_fac,B,k3);

    interp2dVector(&B[0],x0+k3[0],y0+k3[1],z0+k3[2],
            nR_Bfield,nZ_Bfield,BfieldGridR,BfieldGridZ,
            BfieldR,BfieldZ,BfieldT);
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
