
#include <cassert>

template< typename T >
class tensor
{
  public:

  //tensor( T const *data, std::vector< long long unsigned int > const dims );
  tensor( T const *data, long long unsigned int const *dims, int n_dims );

  T get( long long unsigned int *coordinates );

  void set( T val, long long unsigned int *coordinates );

  long long unsigned int const *get_dims();

  /* Ahoy, Captain!!! Make sure this number isn't too small for the problems!!! */
  static int constexpr n_dims_arbitrary_max = 8;

  const int n_dims;

  T const *data;

  long long unsigned int dims[ n_dims_arbitrary_max ];

  long unsigned int offset_factors[ n_dims_arbitrary_max ];
};

template< typename T >
long long unsigned int const *
tensor< T >::get_dims() { return dims; }

/* leading dimension comes first: zyx order */
template< typename T >
tensor< T >::tensor( T const *data, long long unsigned int const *dims_init, int n_dims )
  :
  data( data ),
  n_dims( n_dims )
{ 
    assert( n_dims <= n_dims_arbitrary_max );

    /* populate dims */
    for( int i = 0; i < n_dims; i++ ) dims[ i ] = dims_init[ i ];

    offset_factors[ n_dims - 1 ] = 1;

    for( int i = n_dims - 1; i > 0; i-- )
    {
      offset_factors[ i - 1 ] = offset_factors[ i ] * dims[ i ];
    }
}

/* leading dimension comes first, zyx access */
template< typename T >
T tensor< T >::get( long long unsigned int *coordinates )
{
  int offset = 0;

  /* skip the final element of coordinates  */
  for( int i = 0; i < n_dims; i++ )
  {
    offset += coordinates[ i ] * offset_factors[ i ];
  }

  return data[ offset ];
}

/* Captain! This one has not been tested */
/* leading dimension comes first, zyx access */
template< typename T >
void tensor< T >::set( T val, long long unsigned int *coordinates )
{
  int offset = 0;

  for( int i = 0; i < n_dims; i++ )
  {
    offset += coordinates[ i ] * offset_factors[ i ];
  }

  data[ offset ] = val;
}






