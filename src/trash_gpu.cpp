#include <iostream>
#include <thrust/execution_policy.h>
#include "trash_gpu.h"

operate_boris::operate_boris( )
{
}

CUDA_CALLABLE_MEMBER    
void operate_boris::operator()( std::size_t index )
{
  data[ index ] = data[ index ] * 2;
}

run_boris::run_boris()
  :
  operate(),
  array( 5 * 1e8, 0 )
{
  for( int i = 0; i < array.size(); i++ )
  {
    array[ i ] = i;
  }

  operate.data = array.data();
}

void run_boris::run()
{
  thrust::counting_iterator<std::size_t> particle_iterator_start( 0 );

  thrust::counting_iterator<std::size_t> particle_iterator_end( array.size() );

  thrust::for_each( thrust::device,
      particle_iterator_start,
      particle_iterator_end,
      operate );

  /*
  std::cout << "Ahoy Captain!" << std::endl;
  for( int i = 0; i < array.size(); i++ )
  {
    std::cout << array[ i ] << std::endl;
  }
  */
}
