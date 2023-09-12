#include "random_uniform_numbers.h"

__device__
void
random_uniform_numbers::device_rand_init( long n_particles )
{

  thrust::counting_iterator< std::size_t > start( 0 );
  thrust::counting_iterator< std::size_t > end( n_particles );
  thrust::for_each( thrust::device, 
                    start,
                    end,
                    curandInitialize<>( &particle_state.front(), 0 ) );
}

__host__ __device__
random_uniform_numbers::random_uniform_numbers( long n_particles )
  :
  particle_state( n_particles )
{
  #if USE_CUDA

  device_rand_init( n_particles );
                      
  #else

  std::random_device randDeviceInit;

  for( long i = 0; i < n_particles; ++i )
  {
    std::mt19937 s0(randDeviceInit());
    particle_state[i] = s0;
  }

  #endif
}

float random_uniform_numbers::operator()( int particle_index )
{
  #if USE_CUDA

  return curand_uniform( &particle_state[ particle_index ] );

  #else

  return std::uniform_real_distribution< float >( 0.0, 1.0 )( particle_state[ particle_index ] );

  #endif
}
