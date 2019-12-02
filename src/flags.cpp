#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "flags.hpp"

//CUDA_CALLABLE_MEMBER
//Flags::Flags() : USE_IONIZATION{0} {};
void Flags::initialize_flags(libconfig::Config &cfg) 
{
  USE_IONIZATION = getVariable_cfg<int> (cfg,"flags.USE_IONIZATION");
}
