#include "catch2/catch_all.hpp"
#include "slow_math.h"
#include "constants.h"

TEST_CASE( "slow math tests" )
{
  SECTION( "rotation_matrix" )
  {
    double angle = gitr_constants::pi / 3;

    std::vector< double > rotation_matrix = generate_3d_rotation_matrix( angle );

    REQUIRE( rotation_matrix[ 0 ] == std::cos( angle ) );
    REQUIRE( rotation_matrix[ 1 ] == 0 );
    REQUIRE( rotation_matrix[ 2 ] == std::sin( angle ) );
    REQUIRE( rotation_matrix[ 3 ] == 0 );
    REQUIRE( rotation_matrix[ 4 ] == 1 );
    REQUIRE( rotation_matrix[ 5 ] == 0 );
    REQUIRE( rotation_matrix[ 6 ] == -1 * std::sin( angle ) );
    REQUIRE( rotation_matrix[ 7 ] == 0 );
    REQUIRE( rotation_matrix[ 8 ] == std::cos( angle ) );
  }

  /* matrix dot vector */
  SECTION( "dot_0" )
  {
    std::vector< double > matrix =
    {
      12, 20, 3, 4, 5,
      8, 3, 5, 2, 9,
      6, 23, 29, 2, 6,
      11, 17, 19, 5, 13
    };

    std::vector< double > vector =
    { 4, 7, 2, 3, 5 };

    std::vector< double > gold =
    { 231, 114, 279, 281 };

    int rows = 4;
    int cols = 5;

    std::vector< double > test = slow_dot_product( matrix, vector, rows, cols );

    REQUIRE( test == gold );
  }

  /* vector dot vector */
  SECTION( "dot_1" )
  {
    std::vector< double > matrix =
    { 2, 4, 6, 8, 10 };

    std::vector< double > vector =
    { 9, 7, 5, 3, 11 };

    std::vector< double > gold =
    { 210 };

    int rows = 1;
    int cols = 5;

    std::vector< double > test = slow_dot_product( matrix, vector, rows, cols );

    REQUIRE( test == gold );
  }

  SECTION( "cross" )
  {
    std::vector< double > v0{ 3, 5, 17 };

    std::vector< double > v1{ 11, 8, 15 };

    std::vector< double > gold{ -61, 142, -31 };

    std::vector< double > test = slow_cross_product( v0, v1 );

    REQUIRE( test == gold );
  }
}

TEST_CASE( "rmse comparer tests" )
{
  SECTION( "root_squared" )
  {
    std::vector< double > v0{ 4, 10, 30, 50 };

    std::vector< double > v1( 4 );

    for( int i = 0; i < v0.size(); i++ )
    {
      v1[ i ] = v0[ i ] + 4;
    }

    double error = root_mean_squared_error( v0, v1 );

    REQUIRE( error == 4 );
  }

  SECTION( "rmse_based" )
  {
    int size = 101;

    std::vector< double > v0( size );

    std::vector< double > v1( size );

    double tolerance = 1e-4;

    double element_difference = 0.01;

    for( int i = 0; i < v0.size(); i++ ) v0[ i ] = i;

    for( int i = 0; i < v0.size(); i++ ) v0[ i ] = v1[ i ] + element_difference;

    /* total root error is the sqrt( length * difference ) */
    REQUIRE( rmse_based_comparison( v0, v1, tolerance ) != true );
  }
}
