#include "slow_math.h"

/* calculate root squared error of vectors */
double root_mean_squared_error( std::vector< double > const &v0, 
                           std::vector< double > const &v1 )
{
  assert( v0.size() == v1.size() );

  double squared_sum = 0;

  for( int i = 0; i < v0.size(); i++ )
  {
    double diff = v0[ i ] - v1[ i ];
    squared_sum += diff * diff; 
  }

  return std::sqrt( squared_sum / v0.size() );
}

/* this function compares vectors for approximate equality */
bool rmse_based_comparison( std::vector< double > const &v0, 
                      std::vector< double > const &v1,
                      double const tolerance )
{
  /* compare absolute value */
  auto abs_compare =
  []( double const a, double const b ) -> bool
  { return ( std::abs( a ) < std::abs( b ) ); };

  double max_element =
  std::max( std::abs( *std::max_element( v0.begin(), v0.end(), abs_compare ) ),
            std::abs( *std::max_element( v1.begin(), v1.end(), abs_compare ) ) );

  double normalizing_scale_factor = std::max( 1.0, max_element );

  double normed_rse = root_mean_squared_error( v0, v1 ) / normalizing_scale_factor;

  double scaled_tolerance = tolerance * std::sqrt( v0.size() );

  std::cout << "normed_rse: " << normed_rse
            << " scaled_tolerance: " << scaled_tolerance
            << std::endl;

  return normed_rse < scaled_tolerance;
}

/* call this as a constexpr */
/* generate a row major rotation matrix that rotates "angle" about the y-axis */
std::vector< double > generate_3d_rotation_matrix( double radians_about_y )
{
  /* allocate a 3x3 matrix */
  std::vector< double > const rotation_matrix
  { 
    std::cos( radians_about_y ), 0, std::sin( radians_about_y ),
    0,                           1,                           0,
    -1 * std::sin( radians_about_y ), 0, std::cos( radians_about_y )
  };

  return rotation_matrix;
}

/* dot product function for row major matrix */
std::vector< double > slow_dot_product( std::vector< double > const &matrix,
                                        std::vector< double > const &vector,
                                        int const rows,
                                        int const cols )
{
  assert( matrix.size() == rows * cols );

  std::vector< double > dot_product( rows, 0 );
  
  for( int i = 0; i < rows; i++ )
  {
    for( int j = 0; j < cols; j++ )
    {
      dot_product[ i ] += matrix[ i * cols + j ] * vector[ j ];
    }
  }

  return dot_product;
}

std::vector< double > slow_cross_product( std::vector< double > const &v0,
                                          std::vector< double > const &v1 )
{
  /* cross product is an operation only defined on vectors in 3d Euclidean vector space */
  assert( v0.size() == 3 );
  assert( v1.size() == 3 );

  std::vector< double > crossed( 3 );

  /* 23 - 32 */
  /* 12 - 21 */
  crossed[ 0 ] = v0[ 1 ] * v1[ 2 ] - v0[ 2 ] * v1[ 1 ];

  /* 31 - 13 */
  /* 20 - 02 */
  crossed[ 1 ] = v0[ 2 ] * v1[ 0 ] - v0[ 0 ] * v1[ 2 ];

  /* 12 - 21 */
  /* 01 - 10 */
  crossed[ 2 ] = v0[ 0 ] * v1[ 1 ] - v0[ 1 ] * v1[ 0 ];

  return crossed;
}
