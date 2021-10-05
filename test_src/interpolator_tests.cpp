#include "catch2/catch_all.hpp"
#include "interpolator.h"

/*

s = cell_dimension
d = dimension of interest
fraction = domain and range fraction
low = lower bound bin value
high = upper bound bin value

*/
std::vector< double > generate_basis_cell_1( int s, 
                                             int d,
                                             double fraction,
                                             double low,
                                             double high )
{
  std::vector< double > basis_cell( 1 << s );

  std::cout << "generate_basis_cell_1:" << std::endl;
  std::cout << s << "dimensional-space split across dimension " << d << std::endl;

  for( int i = 0; i < basis_cell.size(); i++ )
  {
    /* if d is 0, use the low. If it is 1, use the high */
    if( ( i >> d ) & 0x1 ) basis_cell[ i ] = high ;

    else basis_cell[ i ] = low;

    std::cout << "basis_cell[ " << i << " ] = " << basis_cell[ i ] << std::endl;
  }

  std::cout << std::endl;

  return basis_cell;
}

/* Ahoy! Build a new basis cell */
std::vector< double > generate_basis_cell_0( int len, int d )
{
  int test = 8;
  int most_significant = ( ( test & ( 0x1 << 3 ) ) >> 3  );
  std::cout << "Ahoy, Captain! most_significant bit: " << most_significant << std::endl;

  std::vector< double > v( 1 << d );

  double len_squared = len * len;

  /* In what order does this interpolate stuff? */
  for( int i = 0; i < 1 << d; i++ )
  {
    double sum = 0;

    for( int j = 0; j < d; j++ )
    {
      sum += ( ( ( i >> j ) & 0x1 ) * len_squared );

      std::cout << ( ( i >> j ) & 0x1 );
    }

    v[ i ] = std::sqrt( sum );

    std::cout << ": " << v[i] << std::endl;
  }

  return v;
}

std::vector< double > embed_basis_cell( std::vector< int > const &d_len,
                                        std::vector< int > const &target_index_tuple,
                                        std::vector< double > const &basis_cell )
{
  /* number of dimensions is d_sizes.size(), 2 ^ n_dims = basis_cell.size(), assert that */
  assert( d_len.size() == target_index_tuple.size() );

  assert( ( 1 << d_len.size() ) == basis_cell.size() );

  /* transform d_len into d_len accumulated */
  std::vector< int > d_multipliers( d_len.size(), 1 );

  /* Captain! Do you need to reverse the coordinates? I think you might... */
  for( int i = 0; i < d_multipliers.size() - 1; i++ )
  {
    d_multipliers[ i + 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  std::cout << "d_multipliers: " << std::endl;
  for( int i = 0; i < d_multipliers.size(); i++ )
  {
    std::cout << " " << d_multipliers[ i ];
  }
  std::cout << std::endl;

  int total_points = d_multipliers.back() * d_len.back();
  std::cout << "total points: " << total_points << std::endl;

  int base_index = 0;

  for( int i = 0; i < d_multipliers.size(); i++ )
  {
    base_index += d_multipliers[ i ] * target_index_tuple[ i ];
  }

  std::cout << "embed base_index: " << base_index << std::endl;
  
  std::vector< double > haystack( total_points, -1 );

  for( int i = 0; i < basis_cell.size(); i++ )
  {
    int flat_index = base_index;
    //std::cout << "i = " << i << std::endl << std::endl;

    for( int j = 0; j < d_multipliers.size(); j++ )
    {
      //std::cout << "j = " << j << " d_multipliers[ j ] = " << d_multipliers[ j ] << std::endl;

      /* Captain! This is still wrong the problem is this needs to be just logical */
      //std::cout << "i & ( 0x1 << j ) = " << ( i & ( 0x1 << j ) ) << std::endl;

      flat_index += ( ( ( i & ( 0x1 << j ) ) >> j ) * d_multipliers[ j ] );

      //std::cout << "flat_index = " << flat_index << std::endl;
    }

    haystack[ flat_index ] = basis_cell[ i ];
    std::cout << "haystack[ " << flat_index << " ] = " 
              << "basis_cell[ " << i << " ]" << std::endl << std::endl;
  }

  return haystack;
}

/* Captain! A row-major, flat data layout is assumed. Note that. */
TEST_CASE( "multi-dimensional interpolation" )
{
  /* see if embed and fetch are truly inverses like you think */
  /*

  Terminology:

  A discrete function is a pair of finite sets.
  The first set:

  x = { x0, x1, x2, ... , x_n }

  contains domain values "x" and the other set is the co-domain values f:

  f = { f(x0), f(x1), f(x2), ... , f(x_n) }

  A continuous function is the same idea but the first set "x" contains all possible domain
  values and the second set "f" contains the co-domain values for every value in "x".

  A linear interpolation is a method to infer a continuous function from a discrete one:

  Given finite "x" and "f" above, linear interpolation attempts to estimate the value of
  f( p ) given:

  p is not contained in set "x"

  there exists x_lower and x_upper in set "x" such that

  x_lower <= p <= x_upper

  This value is inferred by assuming that the function "f" varies linearly and continuously
  in [ x_lower, x_upper ]. Stated more simply, the function on [ x_lower, x_upper ] is assumed to
  be a line, and any value "p" as defined above will lie on that line somewhere.

  The [ x_lower, x_upper ] range containing "p" above will be referred to as the "basis cell"
  containing "p". For "p" in a 1-dimensional vector space, the basis cell contains 2^1 elements
  that define what "p" is "between".

  according to algebra, "p" occupies a position that is

  ( ( p - x_lower ) / ( x_upper - x_lower ) ) * ( x_upper - x_lower ) distance away from x_lower

  and

  ( ( x_upper - p ) / ( x_upper - x_lower ) ) * ( x_upper - x_lower ) distance away from x_upper

  Linear interpolation assumes 

  p = ( ( p - x_lower ) / ( x_upper - x_lower ) ) * ( x_upper - x_lower ) + x_lower 

  and 

  p = x_upper - ( ( x_upper - p ) / ( x_upper - x_lower ) ) * ( x_upper - x_lower )

  Captain! Here
  


  In general, any point "q" in an n-dimensional vector space will have a basis cell containing
  2^n elements that are vertices of an n-dimensional cube containing the point "p".

  x_lower - x_lower <= p - x_lower <= x_upper - x_lower

  */

  SECTION( "embed and fetch a basis cell" )
  {
    /* Captain! Next, make this multidimensional */
    int basis_len = 10;

    int dims = 3;

    auto v0 = generate_basis_cell_0( basis_len, dims );

    /* points in each dim - the space is stretched out so the math remains the same
       as if it was a uniform grid space */
       /* these are the cell dimensions of the space */
    std::vector< int > d_len{ 6, 4, 10 };

    /* this should result in a spacing of 2 in the domain for each dimension */
    /* these are the float space x domain values */
    std::vector< double > max_values{ 12, 8, 20 };

    /* coordinates where it goes */
    /* target location in cell indices */
    std::vector< int > target_index_tuple{ 3, 2, 5 };

    /* the interpolation set for this should be 8 points starting with the one above */
    /* target location cell indices element wise multiplied with spacing plus some 
       arbitrary offset */
    /* Captain! Express coordinates in terms of the variables in the comment above */
    //std::vector< double > coordinates{ 6.9, 4.9, 10.9 };
    std::vector< double > coordinates{ 7.8, 5.8, 11.8 };

    std::vector< double > spacing( 3 );

    for( int i = 0; i < spacing.size(); i++ )
    {
      spacing[ i ] = max_values[ i ] / double(d_len[ i ]);
    }

    auto data = embed_basis_cell( d_len,
                                  target_index_tuple,
                                  v0 );

    auto v1 = fetch_basis_cell( data,
                                coordinates, 
                                d_len,
                                max_values );

    REQUIRE( v0 == v1 );

    double uniform_fraction = 0.5;

    /* Ahoy! This is the part of the test that does not work */
    /*
    std::vector< double > test_point( coordinates.size() );

    std::cout << "Captain! End of the test:" << std::endl;

    for( int i = 0; i < test_point.size(); i++ )
    {
      test_point[ i ] = 
      ( double( target_index_tuple[ i ] ) + uniform_fraction ) * spacing[ i ];
      std::cout << "test_point[ " << i << " ] = " << test_point[ i ] << std::endl;
    }

    double correct_value = uniform_fraction * v1.back();

    double interpolated_value = interpolate_basis_cell( v1, test_point, d_len, max_values );

    REQUIRE( correct_value == interpolated_value );
    */
  }

  SECTION( "interpolate basis cell" )
  {
    /* Captain! Declare these parameters */
    double fraction = 0.6;

    double low = 10;

    double high = 20;

    /* the domain spacing value - space in between domain bins */
    double spacing = 8;

    /* number of cells that span dimension d, assuming the start domain value is 0 */
    double d_len = 5;

    /* note: this is the maximum domain value */
    double max_value = d_len * spacing;

    double gold = ( 1 - fraction ) * low + fraction * high;
    
    /* create a domain value that maps to the gold interpolated value */
    /* put the point in the final cell */
    double test_domain_value = ( d_len - 1 ) * spacing + fraction * spacing;

    /* number of co-planar test points - only variation is within the norm of the plane dirction
       so all interpolations of points in this plane should produce the same answer */
    int n_test_points = 5; 

    int initial_dim = 1;

    int final_dim = 6;

    for( int s = initial_dim; s < final_dim; s++ )
    {
      std::cout << "dimension: " << s << std::endl;

      /* make the 1d test template constants into n-dimensional vectors for interpolations */
      /* later, make a set of test points with dimension d's domain value being the test
         but vary the other dimensional coordinates to be anywhere in the cell. Set element
         "d" below to the test_domain_value and vary the others. Set them all the same then
         just set the "d" element to test_domain_value in the loop */
      std::vector< double > const test_point( s, test_domain_value );
      std::vector< int > const d_len_s( s, d_len );
      std::vector< double > const max_value_s( s, max_value );

      /* iterate the bits in a point in s-dimensional space */
      /* assume zyx is the bit representation order */
      for( int d = 0; d < s; d++ )
      {
        std::vector< double > const basis_cell =
        generate_basis_cell_1( s, d, fraction, low, high );


        /*

        each test point will have values of
        ( d_len - 1 ) * spacing + i / n_test_points * spacing 
        for i going from 1 to n_test_points - 1.

        Add this dimension to the loop - you should get the same
        answer. Verify correct behavior, then clean up. Also put in a variadic template test case

        */

        for( int i = 1; i < n_test_points; i++ )
        {
          double irrelevant_value = ( d_len - 1 ) * spacing + i / n_test_points * spacing;
          std::vector< double > test_point( s, irrelevant_value );
          test_point[ d ] = test_domain_value;

          /* Captain! Test a point set here, after getting 1 to work */
          double test_value = 
          interpolate_basis_cell( basis_cell, test_point, d_len_s, max_value_s );

          double check = std::abs( gold - test_value ) / test_value;
          REQUIRE( check < 1e-10 );
        }
      }
    }

    //REQUIRE( false );
  }
}





















