#ifndef _HISTORY_
#define _HISTORY_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include <math.h>
#include <vector>

struct history { 
    Particles *particlesPointer;
    int tt;
    int subSampleFac;
    int nP;
    float *histX;
    float *histY;
    float *histZ;
    float *histvx;
    float *histvy;
    float *histvz;
    float *histcharge;

    history(Particles *_particlesPointer, int _tt,int _subSampleFac, int _nP, float *_histX,float *_histY,float *_histZ,
          float *_histvx,float *_histvy,float *_histvz, float * _histcharge) : 
        particlesPointer(_particlesPointer), tt(_tt),subSampleFac(_subSampleFac), nP(_nP), 
        histX(_histX),histY(_histY),histZ(_histZ),histvx(_histvx),histvy(_histvy),histvz(_histvz), histcharge(_histcharge) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const 
    {
       if (tt % subSampleFac == 0)
       {
          histX[(tt/subSampleFac)*nP + indx] = particlesPointer->xprevious[indx];
          histY[(tt/subSampleFac)*nP + indx] = particlesPointer->yprevious[indx];
          histZ[(tt/subSampleFac)*nP + indx] = particlesPointer->zprevious[indx];
          histvx[(tt/subSampleFac)*nP + indx] = particlesPointer->vx[indx];
          histvy[(tt/subSampleFac)*nP + indx] = particlesPointer->vy[indx];
          histvz[(tt/subSampleFac)*nP + indx] = particlesPointer->vz[indx];
          histcharge[(tt/subSampleFac)*nP + indx] = particlesPointer->charge[indx];
       } 

    }
};

#endif
