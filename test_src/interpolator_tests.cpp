#include "test_utils.hpp"
#include "interpolator.h"

/* build basis cell, adapt to interpolate it */
std::vector< double > basis_cell( int len, int d )
{
  std::vector< double > v( 1 << d );

  double len_squared = len * len;

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

/* The tests should form a circle back to a known answer:

generate_basis_cell with known interpolation characteristics <---> interpolate_basis_cell

*/

/* Captain! A row-major, flat data layout is assumed. Note that. */
TEST_CASE( "n-dimensional interpolation" )
{
  /* see if embed and fetch are truly inverses like you think */
  SECTION( "embed and fetch a basis cell" )
  {
    int basis_len = 10;

    int dims = 3;

    auto v0 = basis_cell( basis_len, dims );

    /* points in each dim */
    std::vector< int > d_len{ 6, 4, 10 };

    /* this should result in a spacing of 2 for each dimension */
    std::vector< double > max_values{ 12, 8, 20 };

    /* coordinates where it goes */
    std::vector< int > target_index_tuple{ 3, 2, 5 };

    /* the interpolation set for this should be 8 points starting with the one above */
    std::vector< double > coordinates{ 6.9, 4.9, 10.9 };

    auto data = embed_basis_cell( d_len,
                                  target_index_tuple,
                                  v0 );

    /* Captain! You now need float coordinates that you know will map to the coordinates above */
    auto v1 = fetch_basis_cell( data,
                                coordinates, 
                                d_len,
                                max_values );

    /* print out the data */

    /* get this to work before you do the interpolations */
    REQUIRE( v0 == v1 );
  }

  /* this is the hypotenuse test, relies on the geometric identity */
  SECTION( "interpolate an nd basis cell" )
  {
  }

  /* for this test, create a 4d dataset showing a time series progression of
     a 3d field. Make sure the interpolation places the value in the cell center. Just
     add to all values in the cube? */
  SECTION( "interpolate a 4d dataset, cell center, time average of two values." )
  {
  }
}

/* once these are all done you can go workout and swim!!! */




























