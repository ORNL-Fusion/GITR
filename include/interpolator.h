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
#include <vector>
#include <cassert>
#include <cmath>

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
double interpolate_hypercube( std::vector< double > const &hypercube,
                               std::vector< double > const &coordinates,
                               std::vector< int > const &d_len,
                               std::vector< double > const &max_range )
{
  assert( coordinates.size() == d_len.size() );

  assert( coordinates.size() == max_range.size() );

  std::vector< double > spacing( max_range.size() );

  for( int i = 0; i < spacing.size(); i++ )
  {
    spacing[ i ] = max_range[ i ] / double(d_len[ i ]);
  }

  /* the "upper" and "lower" fractions for each dimension */
  std::vector< double > normalized_fractions( coordinates.size() * 2 );

  for( int i = 0; i < coordinates.size(); i++ )
  {
    /* high fraction, matches with a 1 bit */
    normalized_fractions[ i * 2 + 1 ] =
    ( coordinates[ i ] - ( std::floor( coordinates[ i ] / spacing[ i ] ) * spacing[ i ] ) ) /
    spacing[ i ];

    /* low fraction, matches with a 0 bit */
    normalized_fractions[ i * 2 ] = 1 - normalized_fractions[ i * 2 + 1 ];
  }

  /* sum the value of each vertex weighted properly... */
  double sum = 0;

  for( int i = 0; i < hypercube.size(); i++ )
  {
    double weight = 1;

    /* iterate the bits in "i" */
    for( int j = 0; j < coordinates.size(); j++ )
    {
      weight *= normalized_fractions[ j * 2 + ( ( i >> j ) & 0x1 ) ];
    }

    sum += weight * hypercube[ i ];
  }

  return sum;
}


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
std::vector< double > fetch_hypercube( std::vector< double > const &data,
                                       std::vector< double > const &coordinates,
                                       std::vector< int > const &d_len,
                                       std::vector< double > const &max_range )
{
  assert( coordinates.size() == d_len.size() );

  assert( coordinates.size() == max_range.size() );

  /* d_multipliers[ i ] indicates the stride of dimension "i" */
  std::vector< int > d_multipliers( d_len.size(), 1 );

  for( int i = 0; i < d_multipliers.size() - 1; i++ )
  {
    d_multipliers[ i + 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  assert( data.size() == d_multipliers.back() * d_len.back() );

  std::vector< double > spacing( max_range.size() );

  for( int i = 0; i < spacing.size(); i++ )
  {
    spacing[ i ] = max_range[ i ] / double(d_len[ i ]);
  }

  /* find the index of the first vertex in the hypercube */
  int corner_vertex_index = 0;

  for( int i = 0; i < coordinates.size(); i++ )
  {
    corner_vertex_index += 
    ( std::floor( coordinates[ i ] / spacing[ i ] ) * d_multipliers[ i ] );
  }

  std::vector< double > hypercube( 1 << d_len.size() );

  /* find the indices of the other vertices in the hypercube by offsetting from the
     first vertex */
  for( int i = 0; i < hypercube.size(); i++ )
  {
    int flat_index = corner_vertex_index;

    for( int j = 0; j < coordinates.size(); j++ )
    {
      flat_index += ( ( ( i >> j ) & 0x1 ) * d_multipliers[ j ] );
    }

    hypercube[ i ] = data[ flat_index ];
  }

  return hypercube;
}

/* class to encapsulate interpolation functionality */
class interpolated_field
{
  public:

    interpolated_field( std::vector< double > const field,
                        std::vector< int > const d_len,
                        std::vector< double > const max_range )
      :
      field( field ),
      d_len( d_len ),
      max_range( max_range )
      { }


    double operator()( std::vector< double > const coordinates );

  private:

    std::vector< double > const field;

    std::vector< int > const d_len;

    std::vector< double > const max_range;
};

/* Captain! These do not make sense to be passed by value */
double interpolated_field::operator()( std::vector< double > const coordinates )
{
  std::vector< double > hypercube = fetch_hypercube( field, coordinates, d_len, max_range );

  double interpolated_value = interpolate_hypercube( hypercube, coordinates, d_len, max_range );

  return interpolated_value;
}
