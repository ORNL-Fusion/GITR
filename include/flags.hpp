#ifndef _FLAGS_
#define _FLAGS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "libconfig.h++"
#include "utils.h"

class Flags : public ManagedAllocation
{
  private:
  public:
   int USE_IONIZATION;
   CUDA_CALLABLE_MEMBER
   Flags() : USE_IONIZATION{0} {};
   void initialize_flags(libconfig::Config &cfg);
};
#endif
