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
   const bool USE_IONIZATION;
   const bool FIXED_SEEDS;
   const bool USE_ADAPTIVE_DT;
   CUDA_CALLABLE_MEMBER
   Flags(libconfig::Config &cfg) : 
     USE_IONIZATION{initialize_flags(cfg,"USE_IONIZATION")},
       FIXED_SEEDS{initialize_flags(cfg,"FIXED_SEEDS")},
       USE_ADAPTIVE_DT{initialize_flags(cfg,"USE_ADAPTIVE_DT")} {};
   bool initialize_flags(libconfig::Config &cfg, std::string s);
};
#endif
