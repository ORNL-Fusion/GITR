#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "flags.hpp"

//CUDA_CALLABLE_MEMBER
//Flags::Flags() : USE_IONIZATION{0} {};
bool Flags::initialize_flags(libconfig::Config &cfg,std::string s, int default_value) 
{
  std::string base = "flags.";
  //int flag = getVariable_cfg<int> (cfg, base+s);
  int flag;
  if(cfg.lookupValue(base+s, flag))
    {
      std::cout << base+s << " = " << flag << std::endl;
    }
  else
    {
      flag = default_value;
      std::cout << "WARNING: Failed importing " << base+s << ", defaulting value to "<< flag <<  std:: endl;
    }
  if(flag > 0) return true;
  else return false;
}
