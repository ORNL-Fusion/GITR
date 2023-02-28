/* Captain! You must remove vector as a dependency if this is to run on gpus */
#include <vector>

template< typename T >
class tensor
{
  public:

  tensor( T const *data, std::vector< long long unsigned int > const dims );

  T get( std::vector< long long unsigned int > coordinates );

  void set( T val, std::vector< long long unsigned int > coordinates );

  std::vector< long long unsigned int > const &get_dims();

  T const *data;

  std::vector< long long unsigned int > const dims;

  std::vector< long unsigned int > offset_factors;
};

template< typename T >
std::vector< long long unsigned int > const &
tensor< T >::get_dims() { return dims; }

/* leading dimension comes first: zyx order */
template< typename T >
tensor< T >::tensor( T const *data, std::vector< long long unsigned int > const dims )
  :
  data( data ),
  dims( dims )
{ 
    offset_factors.resize( dims.size() );

    offset_factors.back() = 1;

    for( int i = dims.size() - 1; i > 0; i-- )
    {
      offset_factors[ i - 1 ] = offset_factors[ i ] * dims[ i ];
    }
}

/* leading dimension comes first, zyx access */
template< typename T >
T tensor< T >::get( std::vector< long long unsigned int > coordinates )
{
  assert( coordinates.size() == dims.size() );

  int offset = 0;

  /* skip the final element of coordinates  */
  for( int i = 0; i < coordinates.size(); i++ )
  {
    offset += coordinates[ i ] * offset_factors[ i ];
  }

  return data[ offset ];
}

/* Captain! This one has not been tested */
/* leading dimension comes first, zyx access */
template< typename T >
void tensor< T >::set( T val, std::vector< long long unsigned int > coordinates )
{
  assert( coordinates.size() == dims.size() );

  int offset = 0;

  for( int i = 0; i < coordinates.size(); i++ )
  {
    offset += coordinates[ i ] * offset_factors[ i ];
  }

  data[ offset ] = val;
}






