#include "catch2/catch_all.hpp"
#include "interpolator.h"

/*

Purpose:

  Create an n-dimensional Euclidean, unit hypercube.

Parameters:

  n: total number of dimensions
  d: index of the only non-uniform dimension
  low: function value at half of the vertices
  high: function value at the other half of the vertices
  returns: n-dimensional hypercube with all points in a plane normal to "d" having the same value

Concrete Example for 3D:

  In 3D, this is a cube of side length 1, defined by 8 vertices in a linear array
  of size 8, in which the index "i" of each vertex in the array is represented by 3 bits,
  in order of most significant to least,
  representing a Euclidean coordinate ( z, y, x ).

Intent:

  parameter "d" indicates which dimension (x, y, z, ..., n?) is the "varying dimension".
  An n-dimensional unit hypercube will have 2^n vertices. If they are ordered as described in the
  interpolator.h documentation,
  half of them will reside on the hyperplane normal to dimension d where d = 0 and the other
  half will reside on the hyperplane normal to dimension d where d = 1.

  The interpolated value of any coordinate enclosed in this hypercube is thus equivalent to the
  the value of an equivalent 1-dimensional interpolation only in the varying dimension specified.

*/
std::vector< double > generate_hypercube( int n, 
                                             int d,
                                             double low,
                                             double high )
{
  std::vector< double > hypercube( 1 << n );

  for( int i = 0; i < hypercube.size(); i++ )
  {
    /* if d is 0, use the low. If it is 1, use the high */
    if( ( i >> d ) & 0x1 ) hypercube[ i ] = high ;

    else hypercube[ i ] = low;
  }

  return hypercube;
}

/*

Purpose:

  "Hide" a hypercube at lattice coordinates "hypercube_lattice_coordinates"

Parameters:

  d_len:

    for all "i" in _len, d_len[ i ] specifies the number of hypercubes spanning dimension "i"

  hypercube_lattice_coordinates:

    zero-indexed tuple indicated the corner point in lattice coordinates
    at which to place the hypercube

  hypercube:

    2^n values describing a hypercube with corner point at hypercube_lattice_coordinates in
    lattice coordinates

  returns:

    a flattened data array containing values for n-dimensional continuous function "f"
    at points in an ordered n-dimensional lattice representing a regularly sampled subset
    of the domain. The function "f" equals -1 at each lattice point except the points
    defining the n-dimensional hypercube with the smallest point located at
    hypercube_lattice_coordinates

Intent:

    Hide a specific hypercube in a dataset to test that you can correctly retrieve it

*/
std::vector< double > embed_hypercube( std::vector< int > const &d_len,
                                        std::vector< int > const &hypercube_lattice_coordinates,
                                        std::vector< double > const &hypercube )
{
  /* number of dimensions is d_sizes.size(), 2 ^ n_dims = hypercube.size(), assert that */
  assert( d_len.size() == hypercube_lattice_coordinates.size() );

  assert( ( 1 << d_len.size() ) == hypercube.size() );

  /* d_multipliers[ i ] is the stride in dimension "i" */
  std::vector< int > d_multipliers( d_len.size(), 1 );

  for( int i = 0; i < d_multipliers.size() - 1; i++ )
  {
    d_multipliers[ i + 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  int total_points = d_multipliers.back() * d_len.back();

  int base_index = 0;

  for( int i = 0; i < d_multipliers.size(); i++ )
  {
    base_index += d_multipliers[ i ] * hypercube_lattice_coordinates[ i ];
  }

  std::vector< double > haystack( total_points, -1 );

  for( int i = 0; i < hypercube.size(); i++ )
  {
    int flat_index = base_index;

    for( int j = 0; j < d_multipliers.size(); j++ )
    {
      flat_index += ( ( ( i >> j ) & 0x1 ) * d_multipliers[ j ] );
    }

    haystack[ flat_index ] = hypercube[ i ];
  }

  return haystack;
}

TEST_CASE( "multi-dimensional interpolation" )
{

  /* Given an arbitrary coordinate in the domain of "f", fetch the correct enclosing
     hypercube */
  SECTION( "embed and fetch a basis cell" )
  {
    /* number of lattice divisions in each dimension */
    std::vector< int > d_len{ 6, 4, 10 };

    /* difference between initial and final value in each dimension */
    std::vector< double > max_values{ 12, 8, 20 };

    /* lattice coordinates where a hypercube will be hidden to test retrieval */
    std::vector< int > hypercube_lattice_coordinates{ 3, 2, 5 };

    int dimensions = d_len.size();

    assert( dimensions == d_len.size() );
    assert( dimensions == max_values.size() );
    assert( dimensions == hypercube_lattice_coordinates.size() );

    /* just a random point enclosed in the hypercube */
    std::vector< double > enclosed( dimensions );

    for( int i = 0; i < dimensions; i++ )
    {
      enclosed[ i ] = max_values[ i ] / static_cast< double >( d_len[ i ] )
                      * ( static_cast< double >( hypercube_lattice_coordinates[ i ] ) + 0.5 );
    }

    /* create 3 dimensional hypercube symmetric about the first dimension */
    int dimension_of_symmetry = 1;

    double low = 10;

    double high = 20;

    auto hypercube_in = generate_hypercube( dimensions,
                                  dimension_of_symmetry,
                                  low,
                                  high );

    auto lattice = embed_hypercube( d_len,
                                    hypercube_lattice_coordinates,
                                    hypercube_in );

    auto hypercube_out = fetch_hypercube( lattice,
                                          enclosed, 
                                          d_len,
                                          max_values );

    REQUIRE( hypercube_in == hypercube_out );
  }

  /* Given an n-dimensional hypercube and a point enclosed within, interpolate the value of the
     point */
  SECTION( "interpolate hypercube" )
  {
    /* value for the "upper fraction" described in interpolator.h documentation */
    double fraction = 0.6;

    double low = 10;

    double high = 20;

    /* the domain spacing value - space in between hypercube vertices */
    double spacing = 8;

    /* number of hypercubes that span dimension d, assuming the start domain value is 0 */
    double d_len = 5;

    double max_value = d_len * spacing;
    
    /* create a domain value that maps to the gold interpolated value */
    /* put the point in the final cell */
    /* in a one-dimensional interpolation, */

    /* correct 1d interpolation value of a point located "fraction" of the distance between
       "low" and "high" */
    double gold = ( 1 - fraction ) * low + fraction * high;

    double test_domain_value = fraction * spacing;

    /*
    The following loop interpolates "n_test_points" number of points enclosed in an n-dimensional
    hypercube for "n" in a range. Each test point will interpolate to the same known value.
    */
    int n_test_points = 5; 

    int initial_dim = 1;

    int final_dim = 6;

    /* in an s-dimensional space... */
    for( int s = initial_dim; s < final_dim; s++ )
    {
      std::vector< int > const d_len_s( s, d_len );

      std::vector< double > const max_value_s( s, max_value );

      /* ...create an s-dimensional hypercube... */
      for( int d = 0; d < s; d++ )
      {
        /* ...that is symmetric about dimension "d"... */
        std::vector< double > const hypercube =
        generate_hypercube( s, d, low, high );

        /* ...in this hypercube, interpolate the value of n_test_points number of points and
           check against the known correct value */
        for( int i = 1; i < n_test_points; i++ )
        {
          /* 

          points laying in a plane normal to dimension "d" will all interpolate to the same
          value because the coplanar vertices in the planes bounding "d" have the same values:

          these hypercubes are n-dimensional analogs of:

          low        low + fraction * (low - high)              high
           *-----------------*------------------------------------*

          */
          double irrelevant_value =
          static_cast< double >( i ) / static_cast< double >( n_test_points ) * spacing;

          std::vector< double > test_point( s, irrelevant_value );

          test_point[ d ] = test_domain_value;

          double test_value = 
          interpolate_hypercube( hypercube, test_point, d_len_s, max_value_s );

          double check = std::abs( gold - test_value ) / test_value;

          REQUIRE( check < 1e-15 );
        }
      }
    }
  }

  SECTION( "within a class..." )
  {
    /* value for the "upper fraction" described in interpolator.h documentation */
    double fraction = 0.6;

    double low = 10;

    double high = 20;

    /* the domain spacing value - space in between hypercube vertices */
    double spacing = 8;

    /* number of hypercubes that span dimension d, assuming the start domain value is 0 */
    double d_len = 5;

    double max_value = d_len * spacing;
    
    /* create a domain value that maps to the gold interpolated value */
    /* put the point in the final cell */
    /* in a one-dimensional interpolation, */

    /* correct 1d interpolation value of a point located "fraction" of the distance between
       "low" and "high" */
    double gold = ( 1 - fraction ) * low + fraction * high;

    double test_domain_value = fraction * spacing;

    /*
    The following loop interpolates "n_test_points" number of points enclosed in an n-dimensional
    hypercube for "n" in a range. Each test point will interpolate to the same known value.
    */
    int n_test_points = 5; 

    int initial_dim = 1;

    int final_dim = 6;

    /* in an s-dimensional space... */
    for( int s = initial_dim; s < final_dim; s++ )
    {
      std::vector< int > const d_len_s( s, d_len );

      std::vector< double > const max_value_s( s, max_value );

      /* ...create an s-dimensional hypercube... */
      for( int d = 0; d < s; d++ )
      {
        /* ...that is symmetric about dimension "d"... */
        std::vector< double > const hypercube =
        generate_hypercube( s, d, low, high );

        /* ...in this hypercube, interpolate the value of n_test_points number of points and
           check against the known correct value */
        for( int i = 1; i < n_test_points; i++ )
        {
          /* 

          points laying in a plane normal to dimension "d" will all interpolate to the same
          value because the coplanar vertices in the planes bounding "d" have the same values:

          these hypercubes are n-dimensional analogs of:

          low        low + fraction * (low - high)              high
           *-----------------*------------------------------------*

          */
          double irrelevant_value =
          static_cast< double >( i ) / static_cast< double >( n_test_points ) * spacing;

          std::vector< double > test_point( s, irrelevant_value );

          test_point[ d ] = test_domain_value;

          /* embed the hypercube into the middle of the latice somewhere */
          std::vector< int > hypercube_lattice_coordinates = d_len_s;

          for( int j = 0; j < d_len_s.size(); j++ )
          {
            hypercube_lattice_coordinates[ j ] /= 2;
          }

          auto lattice = embed_hypercube( d_len_s,
                                          hypercube_lattice_coordinates,
                                          hypercube );

          class interpolated_field interpolated_field( lattice, d_len_s, max_value_s );

          double test_value_0 = interpolated_field( test_point );

          double check = std::abs( gold - test_value_0 ) / test_value_0;

          REQUIRE( check < 1e-15 );
        }
      }
    }
  }
}
