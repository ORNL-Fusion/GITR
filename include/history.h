#ifndef _HISTORY_
#define _HISTORY_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include "Particles.h"
#include "Boundary.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct history { 
    Particles *particlesPointer;
    int* tt;
    int nT;
    int subSampleFac;
    int nP;
    gitr_precision *histX;
    gitr_precision *histY;
    gitr_precision *histZ;
    gitr_precision *histv;
    gitr_precision *histvx;
    gitr_precision *histvy;
    gitr_precision *histvz;
    gitr_precision *histcharge;
    gitr_precision *histweight;

    history(Particles *_particlesPointer, int* _tt, int _nT,int _subSampleFac, int _nP, gitr_precision *_histX,gitr_precision *_histY,gitr_precision *_histZ,
          gitr_precision *_histv,gitr_precision *_histvx,gitr_precision *_histvy,gitr_precision *_histvz, gitr_precision * _histcharge, gitr_precision * _histweight) : 
        particlesPointer(_particlesPointer), tt(_tt),nT(_nT),subSampleFac(_subSampleFac), nP(_nP), 
        histX(_histX),histY(_histY),histZ(_histZ),histv(_histv),histvx(_histvx),histvy(_histvy),histvz(_histvz), histcharge(_histcharge), histweight(_histweight) {}

    CUDA_CALLABLE_MEMBER_DEVICE    
void operator()(std::size_t indx) const 
    {  
       int tt0=particlesPointer->tt[indx];
       particlesPointer->tt[indx] = particlesPointer->tt[indx]+1;
       //std::cout << "tt subsamplefac indx, nT " << tt0 << " "<< subSampleFac << " " << indx << " " << nT << std::endl;
       if (tt0 % subSampleFac == 0)
       {
       int indexP = particlesPointer->index[indx];
        int histInd = indexP*(nT/subSampleFac+1) + tt0/subSampleFac;
       //std::cout << "histInd " << histInd << std::endl;
        if(histInd <= (nP*(nT/subSampleFac+1)) && histInd >= 0 && indexP < nP)
        {
	  //std::cout << "inside history " << indx << " " << histInd << " " << particlesPointer->zprevious[indx] << std::endl;
          histX[histInd] = particlesPointer->xprevious[indexP];
          histY[histInd] = particlesPointer->yprevious[indexP];
          histZ[histInd] = particlesPointer->zprevious[indexP];
          histv[histInd] = particlesPointer->v[indexP];
          histvx[histInd] = particlesPointer->vx[indexP];
          histvy[histInd] = particlesPointer->vy[indexP];
          histvz[histInd] = particlesPointer->vz[indexP];
          histcharge[histInd] = particlesPointer->charge[indexP];
          histweight[histInd] = particlesPointer->weight[indexP];
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
