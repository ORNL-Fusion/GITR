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
   const bool USE_PARTICLE_DIAGNOSTICS;
   const bool USE_SHEATH_DENSITY;
   CUDA_CALLABLE_MEMBER
   Flags(libconfig::Config &cfg) : 
     USE_IONIZATION{initialize_flags(cfg,"USE_IONIZATION",1)},
       FIXED_SEEDS{initialize_flags(cfg,"FIXED_SEEDS",1)},
       USE_ADAPTIVE_DT{initialize_flags(cfg,"USE_ADAPTIVE_DT",0)},
       USE_PARTICLE_DIAGNOSTICS{initialize_flags(cfg,"USE_PARTICLE_DIAGNOSTICS",0)},
       USE_SHEATH_DENSITY{initialize_flags(cfg,"USE_SHEATH_DENSITY",0)} {};
   bool initialize_flags(libconfig::Config &cfg, std::string s, int default_value);
};
#endif
