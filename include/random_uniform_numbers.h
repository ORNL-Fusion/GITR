/* conditional header includes */

#if USE_CUDA

  #include "curand.h"
  #include "curand_kernel.h"
  #include "thrust/random.h"
  #include "curandInitialize.h"
  using rand_type = curandState;
  #define __hardware_specifier__d  __device__;
  #define __hardware_specifier__h  __host__;

#else

  #include <random>
  #include "thrust/random.h"
  #include <thrust/binary_search.h>
  #include <thrust/execution_policy.h>
  #include <thrust/functional.h>
  #include <thrust/sequence.h>
  #include <thrust/sort.h>
  #include <thrust/transform.h>
  using rand_type = std::mt19937;
  #define __hardware_specifier__d   ;
  #define __hardware_specifier__h   ;

#endif

/* unconditional header includes */
#include "array.h"

/* Declarations */

/* random number generator for each particle */
/* This makes no attempt to utilize MPI - that will come later */
class random_uniform_numbers
{
  public:

    __hardware_specifier__d __hardware_specifier__h
    random_uniform_numbers( long n_particles = 0 );

  #if USE_CUDA
    __hardware_specifier__d
    void device_rand_init( long n_particles );
#endif
    __hardware_specifier__d
    float operator()( int particle_index );

  private:

    sim::Array< rand_type > particle_state;
};
