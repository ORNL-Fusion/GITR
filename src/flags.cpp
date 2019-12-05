#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "flags.hpp"

//CUDA_CALLABLE_MEMBER
//Flags::Flags() : USE_IONIZATION{0} {};
bool Flags::initialize_flags(libconfig::Config &cfg,std::string s) 
{
  std::string base = "flags.";
  int flag = getVariable_cfg<int> (cfg, base+s);
  if(flag > 0) return true;
  else return false;
}
