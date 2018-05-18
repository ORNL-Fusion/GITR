#ifndef _HISTORY_
#define _HISTORY_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"
#include "math.h"
#include <vector>

struct history { 
    Particles *particlesPointer;
    int tt;
    int nT;
    int subSampleFac;
    int nP;
    float *histX;
    float *histY;
    float *histZ;
    float *histvx;
    float *histvy;
    float *histvz;
    float *histcharge;
    float *histweight;

    history(Particles *_particlesPointer, int _tt, int _nT,int _subSampleFac, int _nP, float *_histX,float *_histY,float *_histZ,
          float *_histvx,float *_histvy,float *_histvz, float * _histcharge, float * _histweight) : 
        particlesPointer(_particlesPointer), tt(_tt),nT(_nT),subSampleFac(_subSampleFac), nP(_nP), 
        histX(_histX),histY(_histY),histZ(_histZ),histvx(_histvx),histvy(_histvy),histvz(_histvz), histcharge(_histcharge), histweight(_histweight) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const 
    {
       if (tt % subSampleFac == 0)
       {
        int histInd = indx*nT/subSampleFac + tt/subSampleFac;
        if(histInd < nP*nT/subSampleFac && histInd >= 0 && indx < nP)
        {
          histX[histInd] = particlesPointer->xprevious[indx];
          histY[histInd] = particlesPointer->yprevious[indx];
          histZ[histInd] = particlesPointer->zprevious[indx];
          histvx[histInd] = particlesPointer->vx[indx];
          histvy[histInd] = particlesPointer->vy[indx];
          histvz[histInd] = particlesPointer->vz[indx];
          histcharge[histInd] = particlesPointer->charge[indx];
          histweight[histInd] = particlesPointer->weight[indx];
        }
        //else
        //{
        //    if(particlesPointer->test[indx] == 0.0)
        //    {
        //        particlesPointer->test[indx] = 1.0;
        //        particlesPointer->test0[indx] =     histInd;
        //       particlesPointer->test1[indx] = indx;
        //       particlesPointer->test2[indx]=nT;
        //       particlesPointer->test3[indx]=subSampleFac;
        //      particlesPointer->test4[indx] =tt;
        //    }
        //} 
       }
    }
};

#endif
