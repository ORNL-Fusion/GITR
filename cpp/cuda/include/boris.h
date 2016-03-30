#ifndef _BORIS_
#define _BORIS_

#include "cudaParticle.h"

const float bx=1.0, by=0, bz=0;

__host__ __device__ 
float getE ( float x ) {

    float E = 0;
    return E;
}

struct move_boris { 

        __host__ __device__ 
                void operator()(cudaParticle &p) const { 
                        p.x += 1.0 + getE(p.x) + bx;
                } 
};

#endif
