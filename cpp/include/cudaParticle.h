#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <cstdlib>

class cudaParticle {
	public:
		float x;

		__host__ __device__ cudaParticle();
};
