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
#ifdef __CUDACC__
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#endif

#ifdef __GNUC__ 
#include <stdlib.h>
#endif
#ifdef USE_CUDA
#else
struct ordering {
    bool operator ()(thrust::pair<int, float> const& a,thrust::pair<int, float> const& b) {
            return (b.second) < (a.second);
	        }
};

struct sortParticles { 
    Particles *particles;
    int nP;
    float weightThreshold;
    int& tt;
    int nDtPerApply;
    int* startIndx;
    int* nActiveParticles;
    int rank;
#if __CUDACC__
    curandState *state;
#else
    std::mt19937 *state;
#endif

    sortParticles(Particles *_particles,int _nP, float _weightThreshold,int& _tt,int _nDtPerApply, int* _startIndx,int* _nActiveParticles,int _rank,
#if __CUDACC__
                curandState *_state)
#else
                std::mt19937 *_state)
#endif
    : particles(_particles),nP(_nP), weightThreshold(_weightThreshold),tt(_tt),nDtPerApply(_nDtPerApply), startIndx(_startIndx), nActiveParticles(_nActiveParticles),rank(_rank),state(_state) {}
    
    CUDA_CALLABLE_MEMBER_DEVICE 
    void operator()(std::size_t indx) const 
    { 
    if(tt%nDtPerApply==0 && tt>0)
    {
       int start = startIndx[rank];
       int nPonRank = nActiveParticles[rank];
       int end = start+nPonRank-1;
       sim::Array<float> weight(nPonRank,0.0); 
       sim::Array<thrust::pair<int,float>> pairs(nPonRank);
       //std::cout << "rank tt n start " << rank<<" "<< tt<<" "<< start<< std::endl;
//#ifdef __CUDACC__
	//#else
	//  std::uniform_real_distribution<float> dist(0.0, 1.0);
        //#endif
       for(int i=0;i<nPonRank;i++)
       {
	//#ifdef __CUDACC__
        //weight[i] =curand_uniform(&state[i]);
	//#else
        //weight[i] =dist(state[start+i]);
	//#endif
	//particles->weight[i] = weight[i];
        weight[i] =particles->weight[start+i];
	pairs[i].first = i;
        pairs[i].second = weight[i];
       //std::cout << "pair "  << " " << pairs[i].first << " " << pairs[i].second << std::endl;
       }
       //sim::Array<float> weight0(weight); 
       thrust::sort(thrust::device,pairs.begin(),pairs.end(),ordering());
       for(int i=0;i<nPonRank;i++)
       {
       //std::cout << "pair "  << i<<" " << pairs[i].first << " " << pairs[i].second << std::endl;
       weight[i] = pairs[i].second; 
       }
       sim::Array<float> weightThresholdA(1,weightThreshold);
       sim::Array<int> lowerBoundIndex(1,0);
       //std::cout << "weights " << " " << weightThresholdA[0] << std::endl;
       thrust::upper_bound(weight.begin(), weight.end(),
                           weightThresholdA.begin(),weightThresholdA.end() , 
        		   lowerBoundIndex.begin(),thrust::greater<float>());
       //std::cout << " min index " << lowerBoundIndex[0] << " " << weight[lowerBoundIndex[0]] << std::endl;
       int nUnderThresh=nPonRank-lowerBoundIndex[0];
       std::cout << " nPartivles under thresh " << nUnderThresh << std::endl;
       int nSwapsNeeded=0;
       int goodParticles = nPonRank-nUnderThresh;
       //std::cout << " n good particles " << goodParticles << std::endl;
       for(int i=0;i<nUnderThresh;i++)
       {
        if(pairs[nPonRank-nUnderThresh+i].first < goodParticles)
        {
          nSwapsNeeded=nSwapsNeeded+1;
        }
       }
       //std::cout << " nSwapsNeeded " << nSwapsNeeded << std::endl;
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
       { particles->swapP(swapBad[i]+start,swapGood[i]+start);
       }
       for(int i=0;i<nPonRank;i++)
       {
         //std::cout << " weight0 " << i << " " << particles->weight[start+i]<< " " << particles->index[i] << std::endl;
       }
       nActiveParticles[rank] = goodParticles;
    }
    } 
};
#endif
#endif
