#ifndef _SORT_
#define _SORT_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

#ifdef __GNUC__ 
#include <stdlib.h>
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

struct ordering {
CUDA_CALLABLE_MEMBER_DEVICE
    bool operator ()(thrust::pair<int, gitr_precision> const& a,thrust::pair<int, gitr_precision> const& b) {
            return (b.second) > (a.second);
	        }
};

struct sortParticles { 
    Particles *particles;
    int nP;
    int* tt;
    int nDtPerApply;
    int* nActiveParticles;

    sortParticles(Particles *_particles,int _nP, int* _tt,int _nDtPerApply, int* _nActiveParticles)
    
		    : particles(_particles),nP(_nP), tt(_tt),nDtPerApply(_nDtPerApply), nActiveParticles(_nActiveParticles) {}
    
    CUDA_CALLABLE_MEMBER_HOST 
    void operator()(std::size_t indx) const 
    { 
    if((tt[0]%nDtPerApply==0) & (tt[0]>0))
    {
       int nPonRank = nActiveParticles[0];
       int end = nPonRank-1;
       sim::Array<gitr_precision> hitWall(nPonRank,0); 
       sim::Array<thrust::pair<int,gitr_precision>> pairs(nPonRank);
       
       for(int i=0;i<nPonRank;i++)
       {
        if(particles->weight[i] < 0.0001)
	{
          particles->hitWall[i] = 2.0;
	}
        hitWall[i] =particles->hitWall[i];
	pairs[i].first = i;
        pairs[i].second = hitWall[i];
       }
       thrust::sort(thrust::device,pairs.begin(),pairs.end(),ordering());
       for(int i=0;i<nPonRank;i++)
       {
       //std::cout << "pair "  << i<<" " << pairs[i].first << " " << pairs[i].second << std::endl;
         hitWall[i] = pairs[i].second; 
       }
       sim::Array<gitr_precision> weightThresholdA(1,1.5);
       sim::Array<int> lowerBoundIndex(1,0);
       //std::cout << "weights " << " " << weightThresholdA[0] << std::endl;
       thrust::upper_bound(hitWall.begin(), hitWall.end(),
                           weightThresholdA.begin(),weightThresholdA.end() , 
        		   lowerBoundIndex.begin());
       //std::cout << " min index " << lowerBoundIndex[0] << " " << hitWall[lowerBoundIndex[0]] << std::endl;
       int nUnderThresh=nPonRank-lowerBoundIndex[0];
       //std::cout << " nPartivles under thresh " << nUnderThresh << std::endl;
       int nSwapsNeeded=0;
       int goodParticles = nPonRank-nUnderThresh;
       std::cout << " n good particles " << goodParticles << std::endl;
       for(int i=0;i<nUnderThresh;i++)
       {
        if(pairs[nPonRank-nUnderThresh+i].first < goodParticles)
        {
          nSwapsNeeded=nSwapsNeeded+1;
        }
       }
       //std::cout << " nSwapsNeeded " << nSwapsNeeded << std::endl;
       if(nSwapsNeeded > 0)
       {
         sim::Array<int> swapBad(nSwapsNeeded),swapGood(nSwapsNeeded);
         int ind=0;
         for(int i=0;i<nUnderThresh;i++)
         {
          if(pairs[nPonRank-nUnderThresh+i].first < goodParticles)
          {
            swapBad[ind] = pairs[nPonRank-nUnderThresh+i].first;
         //std::cout << " swapBad " << ind <<" " << swapBad[ind]<<" " << "weight0[swapBad[ind]]" << std::endl;
            ind=ind+1;
          }
         }
         ind=0;
         //std::cout << "swap good going from 0 to " << goodParticles << std::endl;
         for(int i=0;i<goodParticles;i++)
         {
         //std::cout << " swapGood1 " << i <<" "<<ind<< " "  << std::endl;
         //std::cout << " pairs[i].first " << pairs[i].first <<" "<< " "  << std::endl;
          if(pairs[i].first > goodParticles-1)
          {
            swapGood[ind]=pairs[i].first;
         //std::cout << " swapGood " << ind <<" " << swapGood[ind]<<" " << "weight0[swapGood[ind]]" << std::endl;
            ind=ind+1;
          }
         }
         for(int i=0;i<nSwapsNeeded;i++)
         { particles->swapP(swapBad[i],swapGood[i]);
         }
       //for(int i=0;i<nPonRank;i++)
       //{
       //  std::cout << " weight0 " << i << " " << particles->hitWall[i]<< " " << particles->index[i] << std::endl;
       //}
       }
       nActiveParticles[0] = goodParticles;
    }
    } 
};
#endif
