/*

Mathematical Context:

  Terminology:

  A discrete function is a pair of finite sets.
  The first set:

  x = { x0, x1, x2, ... , x_n }

  i < j ---> x[ i ] < x[ j ]

  x[ i ] - x[ i - 1 ] = constant

  contains domain values "x" and the other set is the co-domain values f:

  f = { f(x0), f(x1), f(x2), ... , f(x_n) }

  A continuous function is the same idea but the first set "x" contains all possible domain
  values and the second set "f" contains the co-domain values for every value in "x".

  A linear interpolation is a method to infer a continuous function from a discrete one:

  Given finite "x" and "f" above, linear interpolation attempts to estimate the value of
  f( p ) given:

  p itself is not contained in set "x", but it is in between two elements:

  there exists x_lower and x_upper in set "x" such that

  x_lower <= p <= x_upper

  Assumption:

  For all valid "i" in "x": The Cartesian points

  ( x[ i ], f( x[ i ] ) )

  and

  ( x[ i + 1 ], f( x[ i + 1 ] ) )

  are samples of a continuous function and that function is a straight line on the domain
  interval [ x[ i ], x[ i + 1 ] ] 

  Thus, the value f( p ) of interest must lie along that line:

  f( p ) = a * p + b

  Since this line is within a defined domain interval, "b" can be safely set to 0:

  f( p ) = a * p

  Recall that there exists "x_lower" and "x_upper" in set "x" such that

  x_lower <= p <= x_upper

  where

       ( f( x_upper ) - f( x_lower ) )
  a = --------------------------------- 
            ( x_upper - x_lower )

  but

       ( x_lower * ( x_upper - p ) + x_upper * ( p - x_lower ) ) 
  p = ----------------------------------------------------------- 
                         ( x_upper - x_lower )
  
  so

  
                ( x_lower * ( x_upper - p )             x_upper * ( p - x_lower ) ) 
  f( p ) = a * ------------------------------  +   a * ------------------------------- 
                      ( x_upper - x_lower )                  x_upper - x_lower 

  and from linearity of f:

                               x_upper - p                                p - x_lower
  f( p ) = a * f( x_lower ) * ------------------  +  a * f( x_upper ) *  -------------------
                               x_upper - x_lower                          x_upper - x_lower

  
  Due to equal spacing of the elements in "x", we have:

  there exists an "i" index in "x" such that

  x_lower = x[ i ]
  x_upper = x[ i + 1 ]

  and the matching co-domain values:

  f( x_lower ) = f[ i ]
  f( x_upper ) = f[ i + 1 ]
                               
  due to the spacing of the x[ i ] defined above, x_lower and x_upper can be calculated given "p"
  as follows:

  spacing = x[ j + 1 ] - x[ j ] for any valid index j into "x" 

  x_lower = x[ i ]

  i = std::floor( p / spacing )

  x[ i ] = i * spacing

  thus 

  one-dimensional function f( p, y_const ) for a constant "y_const" is a weighted sum of 
  f( x_lower ) and f( x_upper ) where the
  term applied to f( x_lower ) is called the "lower_fraction"
  and the term applied to
  f( x_upper ) will be referred to as the "upper_fraction"

  f( p, y ) = lower_fraction * f( x_lower, y ) * upper_fraction * f( x_upper, y )

  This can be applied inductively for multidimensional interpolation for calculating for
  input coordinate vector "v":

  v = ( x, y, z, ... )

  f( v ) = f( x = v0, y = v1, z = v2, ... )

  where

  x_lower <= x <= x_upper
  y_lower <= y <= y_upper
  z_lower <= z <= z_upper
  .
  .
  .

  A hypercube in n-dimensions is a set of points defined by all combinations:

  f( x from { x_lower, x_upper },
     y from { y_lower, y_upper },
     z from { z_lower, z_upper }, 
     ... ,
     n from { n_lower, n_upper } )

  This set contains 2^n elements, 

  hypercube = { f(v0), f(v1), ... , f( vn ) }

  each "f(v_)" indexed by a number "i" represented by "n" bits stored
  in "n_bits", where
  a '0' at n_bits[ j ] indicates that dimension "j" in the "v_" value was the "lower" boundary 
  while a '1' bit indicates dimension "j" of the vector "v_" has the "upper" boundary value

  In general, any point "v" in an n-dimensional vector space will have a basis cell containing
  2^n elements that are vertices of an n-dimensional cube containing the point "v".

  An interpolation of f( v ) requires two steps:

  1. Fetching the enclosing hypercube containing "v"
  2. Interpolating the value of f( v )

  f( v ) = sum( w_i * f( b_i ) )

  where 

  w_i is a scalar weight term and
  b_i is the i-th ordered hypercube element

  and "i" can be used to get the value b_i while the "n" bits used to represent "i" can 
  simultaneously be used to encode the "n" values that must be multiplied to get "w_i".

*/
/*

Programmatic Representation:

  "f" are the values of a continuous n-dimensional function sampled at points in the domain.
  This set of points is referred to as the "discrete domain set" and is a proper subset of
  the function's domain, an infinite set.

  The domain set can be spatially represented as the vertices of a finite, n-dimensional lattice.

  For any n-dimensional "coordinate" in the continuous domain, the nearest 2^n lattice points 
  form an
  n-dimensional hypercube. This enclosing hypercube is used to interpolate the function's value
  at "coordinate". The 2^n vertices of this cube can be indexed identically to the 2^n vertices
  of an n-dimensional unit hypercube.

  The n-dimensional unit hypercube has 2^n vertices.

  Each element "i" in the index set:

  { 0, 1, 2, ..., ( 2^n ) - 1 }

  indexes one of the 2^n vertices.

  The binary representation of "i" is a bit string of
  n-bits where the final bit b_1 is the least significant and 
  the first bit b_n is the most significant.
  (Big Endian bit order):

  { bit_n, bit_n-1, ..., bit_z, bit_y, bit_x, ... , bit_2, bit_1 }

  which incidentally can be interpreted as a spacial n-tuple and thus vertex of the hypercube:

  { dim_n, dim_n-1, ..., dim_z, dim_y, dim_x, ... , dim_2, dim_1 }

  A binary n-tuple of this form is associated with each vertex "i" and will be referred to as
  "n-tuple_i".

*/
/*

Final Calculation:

  1. To interpolate a point "p", obtain the function values at each of the 2^n vertices 
     of the enclosing hypercube.

  2. Order the values according to the n-tuple_i value for each value's associated vertex
     in a list called "v" with elements "v_i"

  In the multidimensional case, the interpolation is represented as a sum:

  f( coordinate ) = sum over i: ( w_i * v_i )

  and w_i = 
  product over each bit "b" in n-tuple_i:
  ( lower_fraction if b is '0' or upper_fraction if b is '1' ) 

*/

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "flat_array.h"

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

template< typename T >
class interpolated_field : public tensor< T >
{
  public:

    interpolated_field( T const *data,
                        long long unsigned int const *dims,
                        T const *max_range_init,
                        T const *min_range_init,
                        int n_dims_init )
      :
      tensor< T >( data, dims, n_dims_init )
      { 
        data_size = 1;

        /* Captain! This should be n_bins + 1 */
        for( int i = 0; i < this->n_dims; i++ ) data_size *= dims[ i ];
        for( int i = 0; i < this->n_dims; i++ ) max_range[ i ] = max_range_init[ i ];
        for( int i = 0; i < this->n_dims; i++ ) min_range[ i ] = min_range_init[ i ];

        /* Captain!!! This should be number of bins, not dims! */
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

    /* Captain! Create some private methods here! Move the fetch and interpolate functions
       into here. They should not calculate their own offset factors... */

    /* Captain! Does this handle edge cases for points off of the grid? */
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
      std::cout << "Ahoy!" << std::endl;
      for( int i = 0; i < this->n_dims; i++ )
      {
        int reach = 1 << i;

        int step = reach << 1;

        std::cout << "reducing dim " << i << std::endl;

        for( int j = 0; j < 1 << this->n_dims; j += step )
        {
          std::cout << "summing hypercube vertices " << j << " and " << j + reach << std::endl;
          std::cout << hypercube[ j ] << " " << hypercube[ j + reach ] << std::endl;
          std::cout << "pairing with fractions " << i*2 << " and " << i*2+1 << std::endl;
          std::cout << normalized_fractions[i*2] << " " << normalized_fractions[i*2+1] 
                    << std::endl;

          hypercube[ j ] = 
          ( normalized_fractions[ i * 2 ] * hypercube[ j ] + 
          normalized_fractions[ i * 2 + 1 ] * hypercube[ j + reach ] )
          / spacing[ i ];

          std::cout << "result: " << hypercube[ j ] << std::endl;
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
  T hypercube[ 1 << this->n_dims_arbitrary_max ];

  fetch_hypercube( coordinates, hypercube );

  T interpolated_value = interpolate_hypercube( hypercube, coordinates );

  return interpolated_value;
}
