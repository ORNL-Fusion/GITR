#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif
#include <cassert>
#include <cmath>


template< typename T >
class tensor
{
  public:

  //tensor( T const *data, std::vector< long long unsigned int > const dims );
  CUDA_CALLABLE_MEMBER
  tensor( T const *data, long long unsigned int const *dims, int n_dims );

  CUDA_CALLABLE_MEMBER
  T get( long long unsigned int *coordinates );

  CUDA_CALLABLE_MEMBER
  void set( T val, long long unsigned int *coordinates );

  CUDA_CALLABLE_MEMBER
  long long unsigned int const *get_dims();

  /* magic number - bad */
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

/* untested tested */
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

/* Captain! turn this dummy into a marshalling class for the interpolated field */
template< typename T >
class dummy
{
  public:

  // step 1, the constructor needs to take the same arguments as the interpolated field
  CUDA_CALLABLE_MEMBER
  dummy( T const *data_, int n_dims_init )
    :
    data( data_ )
  {}

  // switch out the data pointer
  CUDA_CALLABLE_MEMBER
  void eat( double* data_  )
  { data = data_; }

  double const *data;
};

template< typename T >
class interpolated_field : public tensor< T >
{
  public:

    CUDA_CALLABLE_MEMBER
    interpolated_field( T const *data,
                        long long unsigned int const *dims,
                        T const *max_range_init,
                        T const *min_range_init,
                        int n_dims_init )
      :
      tensor< T >( data, dims, n_dims_init )
      { 
        data_size = 1;

        /* should this be n_bins + 1? */
        for( int i = 0; i < this->n_dims; i++ ) data_size *= dims[ i ];
        for( int i = 0; i < this->n_dims; i++ ) max_range[ i ] = max_range_init[ i ];
        for( int i = 0; i < this->n_dims; i++ ) min_range[ i ] = min_range_init[ i ];

        for( int i = 0; i < this->n_dims; i++ ) 
          spacing[ i ] = ( max_range[ i ] - min_range[ i ] ) / ( T(dims[ i ]) );

      }


    CUDA_CALLABLE_MEMBER
    T operator()( T const *coordinates );

    CUDA_CALLABLE_MEMBER
    void fetch_hypercube( T const *coordinates, T *hypercube );

    CUDA_CALLABLE_MEMBER
    T interpolate_hypercube( T *hypercube,
                             T const *coordinates );

    T max_range[ tensor< T >::n_dims_arbitrary_max ];

    T min_range[ tensor< T >::n_dims_arbitrary_max ];

    T spacing[ tensor< T >::n_dims_arbitrary_max ];

    int data_size;
};

/*

hypercube: n-dimensional hypercube in a flattened array of 2^n vertices

coordinates: domain coordinates for which it is desired to interpolate a function value

d_len:

  for all "i" in d_len, d_len[ i ] specifies the number of equally spaced samples
  spanning domain dimension "i"

max_range:

  for all "i" in max_range, max_range[ i ] specifies the linear extent divided equally
  by d_len[ i ] bins

*/
template< typename T >
T interpolated_field< T >::interpolate_hypercube( T *hypercube,
                                                  T const *coordinates )
{
  /* the "upper" and "lower" fractions for each dimension */
  double normalized_fractions[ this->n_dims_arbitrary_max * 2 ];

  /* break! 4 */
  /* linear iteration over the dimensions in parallel with coordinates. This imples that the
     coordinates here are in zyx order. */
  for( int i = 0; i < this->n_dims; i++ )
  {
    /* high fraction, matches with a 1 bit */
    normalized_fractions[ i * 2 + 1 ] =
    ( coordinates[ i ] - min_range[ i ] 
    - ( std::floor( ( coordinates[ i ] - min_range[ i ]) / spacing[ i ] ) * spacing[ i ] ) )
    ;// / spacing[ i ];

    /* Captain! I think the 1 minus is the thing that is the problem!!! */
    /* low fraction, matches with a 0 bit */
    //normalized_fractions[ i * 2 ] = 1 - normalized_fractions[ i * 2 + 1 ];

    normalized_fractions[ i * 2 ] =
    ( ( ( std::floor( ( coordinates[ i ] - min_range[ i ]) / spacing[ i ] ) + 1 ) * spacing[ i ] )
    - ( coordinates[ i ] - min_range[ i ] ) );// / spacing[ i ];
  }

  /* sum the value of each vertex weighted properly... */
  double sum = 0;

  /* Captain!!! move one level up now */

  /* linear iteration over the hypercube bits in xyz bit mapping and counting from 
     0 to 2^n - 1 */

  /* each xyz cell index represented by the xyz interpretation of the bits in the index itself,
     are iterated upon and consumed linearly in zyx order in parallel to what must be the order
     of the normalized fractions. This implies that normalized_fractions has the same inherent
     ordering in dimensionality as the strides in offset_factors.  */
  /*
  for( int i = 0; i < ( 1 << this->n_dims ); i++ )
  {
    double weight = 1;

    for( int j = 0; j < this->n_dims; j++ )
    {
      weight *= normalized_fractions[ j * 2 + ( ( i >> j ) & 0x1 ) ];
    }

    sum += weight * hypercube[ i ];
  }

  return sum;
  */
      //std::cout << "Ahoy!" << std::endl;
      for( int i = 0; i < this->n_dims; i++ )
      {
        int reach = 1 << i;

        int step = reach << 1;

        //std::cout << "reducing dim " << i << std::endl;

        for( int j = 0; j < 1 << this->n_dims; j += step )
        {
          /*
          std::cout << "summing hypercube vertices " << j << " and " << j + reach << std::endl;
          std::cout << hypercube[ j ] << " " << hypercube[ j + reach ] << std::endl;
          std::cout << "pairing with fractions " << i*2 << " and " << i*2+1 << std::endl;
          std::cout << normalized_fractions[i*2] << " " << normalized_fractions[i*2+1] 
                    << std::endl;
                    */

          hypercube[ j ] = 
          ( normalized_fractions[ i * 2 ] * hypercube[ j ] + 
          normalized_fractions[ i * 2 + 1 ] * hypercube[ j + reach ] )
          / spacing[ i ];

          //std::cout << "result: " << hypercube[ j ] << std::endl;
        }
      }

      return hypercube[ 0 ];
}

/* Coordinates here go from small stride to large stride! */
/* thus, d_len also goes from small stride to large stride offsets. grows like ---> */

/* you simply need to change the calculation of d_len in here and the order of 
   coordinates in the test file */

/*

data: finite lattice in Euclidean n-space of samples of continuous n-dimensional function f 

coordinates: 

  domain coordinates for which it is desired to interpolate a function value
  coordinates are expected in order of dimension with the smallest stride (1) to dimension
  with the largest stride

d_len:

  for all "i" in d_len, d_len[ i ] specifies the number of equally spaced histogram bins
  spanning dimension "i"

max_range:

  for all "i" in max_range, max_range[ i ] specifies the linear extent divided equally
  by d_len[ i ] bins

returns: n-dimensional hypercube in a flattened array of 2^n vertices

*/
template< typename T >
void interpolated_field< T >::fetch_hypercube( T const *coordinates, T *hypercube )
{
  /* find the index of the first vertex in the hypercube */
  /* break! 3 */
  int corner_vertex_index = 0;

  long unsigned int *of = this->offset_factors;

  for( int i = 0; i < this->n_dims; i++ )
  {
    corner_vertex_index += 
    ( std::floor( ( coordinates[ i ] - min_range[ i ] ) / spacing[ i ] ) * this->offset_factors[ i ] );
  }

  /* find the indices of the other vertices in the hypercube by offsetting from the
     first vertex */
  for( int i = 0; i < ( 1 << this->n_dims ); i++ )
  {
    int flat_index = corner_vertex_index;

    for( int j = 0; j < this->n_dims; j++ )
    {
      /* Captain!!! How do your fractions relate to the order of the offset_factors? */

      /* let's examine... */

      /* "i" the linear index of a hypercube cell index contains n_dims bits: 
         these represent xyz order, such that consumption goes zyx in parallel with
         the linear indices of offset_factors */

      /* thus, the order of the hypercube cell coordinates follows 0 ---> 2^n in binary
         where bits represent xyz order: note that the fastest changing dimension in this
         enumeration is the leading dimension! 

         That seems a bit bad for memory access...
         what if we had overlapping memory regions to mitigate these factors? That would be
         awesome... memory cubes where locality is better or square or something... maybe
         nested squares? Could that work? Yeah... that seems like a better data tiling
         layout... 
          */
      flat_index += ( ( ( i >> j ) & 0x1 ) * this->offset_factors[ j ] );
    }

    if( flat_index >= data_size ) assert( flat_index < data_size );

    hypercube[ i ] = this->data[ flat_index ];
  }
}

/* These do not make sense to be passed by value */
template< typename T >
T interpolated_field< T >::operator()( T const *coordinates )
{
  if( data_size == 1 ) return this->data[ 0 ];

  T hypercube[ 1 << this->n_dims_arbitrary_max ];

  fetch_hypercube( coordinates, hypercube );

  T interpolated_value = interpolate_hypercube( hypercube, coordinates );

  return interpolated_value;
}
