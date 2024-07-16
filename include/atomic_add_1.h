#if USE_CUDA >0
#include "cuda_runtime.h"

__device__ double atomicAdd1(double* address, double val);
#endif
