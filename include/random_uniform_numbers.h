/* conditional header includes */

#if USE_CUDA

  #include "curand.h"
  #include "curand_kernel.h"
  #include "thrust/random.h"
  #include "curandInitialize.h"
  using rand_type = curandState;
  #define __hardware_specifier__  __device__;

#else

  #include <random>
  using rand_type = std::mt19937;
  #define __hardware_specifier__  __host__;

#endif

/* unconditional header includes */
#include "array.h"

/* Declarations */

/* random number generator for each particle */
/* This makes no attempt to utilize MPI - that will come later */
class random_uniform_numbers
{
  public:

    __host__ __device__
    random_uniform_numbers( long n_particles = 0 );

    __device__
    void device_rand_init( long n_particles );

    __hardware_specifier__
    float operator()( int particle_index );

  private:

    sim::Array< rand_type > particle_state;
};
