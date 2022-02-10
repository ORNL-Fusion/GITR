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
   const bool USE_CYLSYMM;
   const int USE_PERPDIFFUSION;
   CUDA_CALLABLE_MEMBER
   Flags(libconfig::Config &cfg) : 
     USE_IONIZATION{initialize_flags(cfg,"USE_IONIZATION")},
       FIXED_SEEDS{initialize_flags(cfg,"FIXED_SEEDS")},
       USE_ADAPTIVE_DT{initialize_flags(cfg,"USE_ADAPTIVE_DT")},
       USE_CYLSYMM{initialize_flags(cfg,"USECYLSYMM")},
      USE_PERPDIFFUSION{initialize_int_flags(cfg,"USEPERPDIFFUSION")} {};
   bool initialize_flags(libconfig::Config &cfg, std::string s);
   int initialize_int_flags(libconfig::Config &cfg, std::string s);
};
#endif
