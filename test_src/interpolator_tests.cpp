#include "catch2/catch_all.hpp"
#include "interpolator.h"
//#include "interpolate_3d.h"
#include "interp2d.hpp"
#include "operation.hpp"
#include <iomanip>
#include <chrono>
#include <io.h>

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
                                        std::vector< int > 
                                        const &hypercube_lattice_coordinates,
                                        std::vector< double > const &hypercube )
{
  /* number of dimensions is d_sizes.size(), 2 ^ n_dims = hypercube.size(), assert that */
  assert( d_len.size() == hypercube_lattice_coordinates.size() );

  assert( ( 1 << d_len.size() ) == hypercube.size() );

  /* d_multipliers[ i ] is the stride in dimension "i" */
  std::vector< int > d_multipliers( d_len.size(), 1 );
  /*

  for( int i = 0; i < d_multipliers.size() - 1; i++ )
  {
    d_multipliers[ i + 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  int total_points = d_multipliers.back() * d_len.back();
  */

  for( int i = d_multipliers.size() - 1; i > 0; i-- )
  {
    d_multipliers[ i - 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  int total_points = d_multipliers.front() * d_len.front();

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
  SECTION( "t0" )
  {
    /* number of lattice divisions in each dimension is n_points - 1 */
    //std::vector< int > d_len{ 6, 4, 10 };
    std::vector< int > d_len{ 7, 5, 11 };

    /* difference between initial and final value in each dimension */
    std::vector< double > max_values{ 12, 8, 20 };
    std::vector< double > min_values{ 0, 0, 0 };

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

    std::vector< double > lattice = embed_hypercube( d_len,
                                    hypercube_lattice_coordinates,
                                    hypercube_in );

    class interpolated_field< double >
    field( lattice.data(), d_len.data(), max_values.data(), min_values.data(), dimensions );
    /* Create the hypercube. What will own the pointers? Once this is done, you
       have all the pieces I think. Then you just need to prepare GITR for compiling and
       write out a python file that can convert between the two file types. Then make
       sure that you have working containers */
    std::vector< double > hypercube_out( 1 << dimensions );

    field.fetch_hypercube( enclosed.data(), hypercube_out.data() );

    REQUIRE( hypercube_in == hypercube_out );
  }

  /* Given an n-dimensional hypercube and a point enclosed within, interpolate the value of the
     point */
  SECTION( "t1" )
  {
    /* value for the "upper fraction" described in interpolator.h documentation */
    double fraction = 0.6;

    /* These are f( x ) values */
    double low = 10;

    double high = 20;

    /* number of hypercubes that span dimension d, assuming the start domain value is 0 */
    double d_len = 1;

    double max_value = 1;
    double min_value = 0;
    double spacing = 1;
    
    /* create a domain value that maps to the gold interpolated value */
    /* put the point in the final cell */
    /* in a one-dimensional interpolation, */

    /* correct 1d interpolation value of a point located "fraction" of the distance between
       "low" and "high" */
    double gold = ( 1 - fraction ) * low + fraction * high;

    double test_domain_value = fraction;

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
      // add +1 to d_len because n_points defining your bin boundaries means n_points - 1 bins
      std::vector< int > const d_len_s( s, d_len + 1 );

      std::vector< double > const max_value_s( s, max_value );
      std::vector< double > const min_value_s( s, 0 );

      /* ...create an s-dimensional hypercube... */
      for( int d = 0; d < s; d++ )
      {
        /* ...that is symmetric about dimension "d"... */
        std::vector< double > hypercube =
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

          std::vector< double > hypercube_workspace = hypercube;

          /* note: you are using the hypercube as the entire data lattice just to test
             interpolation functionality independent of hypercube fetching */
          class interpolated_field< double > 
          field( hypercube.data(), d_len_s.data(), max_value_s.data(), min_value_s.data(), s );

          double test_value = 
          field.interpolate_hypercube( hypercube_workspace.data(), test_point.data() );

          double check = std::abs( gold - test_value ) / test_value;

          //std::cout << "gold: " << gold << " test_value: " << test_value << std::endl;

          REQUIRE( std::abs( check ) < 1e-15 );
        }
      }
    }
  }

  SECTION( "t2" )
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

    /*
    The following loop interpolates "n_test_points" number of points enclosed in an n-dimensional
    hypercube for "n" in a range. Each test point will interpolate to the same known value.
    */
    int n_test_points = 5; 

    int initial_dim = 2;

    int final_dim = 9;

    /* in an s-dimensional space... */
    for( int s = initial_dim; s < final_dim; s++ )
    {
      std::vector< int > const d_len_s( s, d_len + 1 );

      std::vector< double > const max_value_s( s, max_value );
      std::vector< double > const min_value_s( s, 0 );

      /* ...create an s-dimensional hypercube... */
      for( int d = 0; d < s; d++ )
      {
        /* ...that is symmetric about dimension "d"... */
        std::vector< double > hypercube =
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
          /* embed the hypercube into the middle of the latice somewhere */
          std::vector< int > hypercube_lattice_coordinates( s, d_len );

          for( int j = 0; j < d_len_s.size(); j++ )
          {
            hypercube_lattice_coordinates[ j ] /= 2;
          }

          double test_point_offset =
          static_cast< double >( i ) / static_cast< double >( n_test_points );

          double irrelevant_value = 
          spacing * ( hypercube_lattice_coordinates[ 0 ] + test_point_offset );

          /* spacing * ( d_len + test_point_offset ) */


          /* new code */
          std::vector< double > enclosed( s, irrelevant_value );

          enclosed[ d ] = max_value_s[ d ] / static_cast< double >( d_len_s[ d ] - 1 )
              * ( static_cast< double >( hypercube_lattice_coordinates[ d ] ) + fraction );
          /* end new code */

          auto lattice = embed_hypercube( d_len_s,
                                          hypercube_lattice_coordinates,
                                          hypercube );

          class interpolated_field< double >
          field( lattice.data(), d_len_s.data(), max_value_s.data(), min_value_s.data(), s );

          double test_value_0 = field( enclosed.data() );

          double check = std::abs( gold - test_value_0 ) / test_value_0;

          REQUIRE( std::abs( check ) < 1e-14 );
        }
      }
    }
  }

  /* test against multidimensional analytical function */
  SECTION( "t3" )
  {
    // params - generate
    double domain_start = -3.5;
    double domain_end = 3.5;

    // the next two are the only base level varying params
    int const n_grid_points = 501;
    int const n_dims = 3;

    auto analytical_function =
    []( std::array< double, n_dims > &real_space ) -> double
    {
      double d = 0;

      for( int j = 0; j < real_space.size(); j++ ) d += real_space[ j ] * real_space[ j ];

      d = 10 * std::cos( std::sqrt( d ) );

      return d;
    };

    // outer loop sweeps over params
    // construct objects based on params
    // inner loop sweeps over all cells

    // dims is n_grid_points in each dimension
    std::unique_ptr< int[] > dims( new int[ n_dims ] );

    for( int i = 0; i < n_dims; i++ ) dims[ i ] = n_grid_points;

    // max and min range for each dimension
    std::unique_ptr< double[] > min_range_init( new double[ n_dims ] );

    for( int i = 0; i < n_dims; i++ ) min_range_init[ i ] = domain_start;

    std::unique_ptr< double[] > max_range_init( new double[ n_dims ] );

    for( int i = 0; i < n_dims; i++ ) max_range_init[ i ] = domain_end;

    // number of cells and spacing between cell boundaries
    int n_cells = n_grid_points - 1;

    double spacing = ( domain_end - domain_start ) / n_cells;

    // offset + cell_corner_point = cell_midpoint
    std::array< double, n_dims > offset = { 0.5 * spacing, 
                                            0.5 * spacing,
                                            0.5 * spacing };

    // begin function: generate_domain_grid - this is only for Tim's implementations

    // allocate grids for Tim's interpolator - the data is strided along dims
    std::unique_ptr< double[] > grid_ptr( new double[ n_dims * n_grid_points ] );

    // populate square grids
    for( int i = 0; i < n_dims; i++ )
    {
      for( int j = 0; j < n_grid_points; j++ )
      {
        grid_ptr[ i * n_grid_points + j ] = spacing * j + domain_start;
      }
    }

    double *grid_raw = grid_ptr.get();

    // end function: generate_domain_grid

    // function begin: generate_lattice --> generate_ifield()

    // calculate needed size
    int lattice_size = 1;

    for( int i = 0; i < n_dims; i++ )
    {
      lattice_size *= n_grid_points;
    }

    // allocate lattice
    int lattice_cells = 1;

    for( int i = 0; i < n_dims; i++ )
    {
      lattice_cells *= ( n_grid_points - 1 );
    }

    std::unique_ptr< double[] > lattice_ptr( new double[ lattice_size ] );

    double *lattice_data = lattice_ptr.get();

    interpolated_field< double > 
    lattice( lattice_data, 
             dims.get(), 
             max_range_init.get(), 
             min_range_init.get(), 
             n_dims );

    // end function generate_lattice()

    // begin function populate_lattice()

    // populate lattice with data
    round_robin_nd< n_dims, n_grid_points > lattice_indexer;

    for( int i = 0; i < lattice_size; i++ )
    {
      /* grid-space */
      lattice_indexer.back_spin();
      auto lattice_gridpoint = lattice_indexer.get_indices();

      /* Captain! calculate indices and check outcome */
      // iterate over the coordinates
      for( int j = 0; j < n_dims; j++ )
      {
      }


      //if( i % 1000  == 0 ) std::cout << i << "/" << lattice_size << std::endl;

      /*
      std::cout << "grid-space" << std::endl;
      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      std::cout << " " << lattice_gridpoint[ j ]; std::cout << std::endl;

      std::cout << "real-space" << std::endl;
      */
      std::array< double, n_dims > real_space;

      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      real_space[ j ] = lattice_gridpoint[ j ]*spacing+domain_start;

      /*
      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      std::cout << " " << real_space[ j ]; std::cout << std::endl;
      */

      /* calculate function value at real-space coordinate */
      double analytical_value = analytical_function( real_space );

      /* store at lattice coordinate "i" */
      lattice.set( analytical_value, lattice_gridpoint.data() );

      double stored_analytical_value = lattice.get( lattice_gridpoint.data() );

      REQUIRE( analytical_value == stored_analytical_value );
    }

    /* function end: populate_lattice */

    // function begin: generate_test_results()

    /* checking begins */

    int legacy_correct = 0;

    int new_correct = 0;

    round_robin_nd< n_dims, n_grid_points - 1 > cell_indexer;

    // columns: cell_index x y z f_analytical f_old f_new
    // cell_index is the "time column"
    // only stack when one variable is fixed, for easier graphing
    std::filesystem::path output_file_name = "interpolator_comparison_3d.csv";
    csv_row_stacker< 6 > data_stack_3d( output_file_name );

    /* interpolate a point in each lattice cell of the domain */
    double total_difference = 0;
    for( int i = 0; i < lattice_cells; i++ )
    {
      /* convert linear cell index into n-dimensional cell coordinate */
      cell_indexer.back_spin();

      auto lattice_gridpoint = cell_indexer.get_indices();

      /* convert cell coordinate into xyz realspace coordinate for bottom corner of cell */
      std::array< double, n_dims > real_space;
      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      real_space[ j ] = lattice_gridpoint[ j ]*spacing+domain_start;

      /* convert realspace coordinate into the center of the cell instead of the corner */
      std::array< double, n_dims > real_space_offset;
      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      {
        real_space_offset[ j ] = real_space[ j ] + offset[ j ];
      }

      /* generate results */
      double analytical_offset_value = analytical_function( real_space_offset );

      double interpolated_value = lattice( real_space_offset.data() );

      double legacy_interpolation =
      interp3d(  
      real_space_offset[ 2 ],
      real_space_offset[ 1 ],
      real_space_offset[ 0 ],
      n_grid_points,
      n_grid_points,
      n_grid_points,
      grid_raw,
      grid_raw + n_grid_points,
      grid_raw + ( n_grid_points * 2 ),
      lattice_data );

      /* process results */
      double difference = 
      std::abs( ( interpolated_value - analytical_offset_value ) / analytical_offset_value );

      total_difference += difference;

      double legacy_difference =
      std::abs( ( legacy_interpolation - analytical_offset_value ) / analytical_offset_value );

      if( legacy_difference >= difference ) ++new_correct;

      else ++legacy_correct;

      /* output results */
      /* take a projection */
      if( real_space[ 1 ] == 0 )
      {
        // columns: cell_index_i x y z f_analytical f_old f_new
        std::array< double, 6 > 
        d{ real_space_offset[ 2 ], real_space_offset[ 1 ], real_space_offset[ 0 ],
           analytical_offset_value, legacy_interpolation, interpolated_value };

        std::span< double, 6 > data_span( d );

        data_stack_3d.stack( i, data_span );
      }

      /* print test code - check them */
      if( i % 1000000  == 0 )
      {
        // integer division will not work...
        std::cout << i << "/" << lattice_cells << " = " 
                  << (double)i / (double)lattice_cells * 100 << std::endl;

        std::cout << "new method difference: " << std::setprecision( 10 )
          << difference
          << std::endl;

        std::cout << "legacy method difference: " << std::setprecision( 10 )
          << legacy_difference
          << std::endl;

        std::cout << "legacy: " << legacy_correct << " new: " << new_correct << std::endl;

        std::cout << std::setprecision( 10 ) << "discrepancy: "
          << ( interpolated_value - legacy_interpolation ) / analytical_offset_value
          << std::endl;
      }
    }

    double mean_difference = total_difference / lattice_cells;

    std::cout << "normalized mean difference: " << mean_difference << std::endl;

    REQUIRE( mean_difference < 0.0003 );
  }

  /* performance compare legacy 3d interp and new one */
  SECTION( "t4" )
  {
    // start and end of domain
    double domain_start = -3.5;
    double domain_end = 3.5;
    int const n_grid_points = 401;

    // defined by leftmost boundary value
    int n_cells = n_grid_points - 1;

    // the above are uniform over n_dims dimensions - start with 3
    int const n_dims = 3;

    std::array< double, n_dims > min_range_init = { domain_start, domain_start, domain_start };

    std::array< double, n_dims > max_range_init = { domain_end, domain_end, domain_end};

    double spacing = ( domain_end - domain_start ) / ( n_grid_points - 1 );

    std::array< double, n_dims > offset = { 0.5 * spacing, 
                                            0.5 * spacing,
                                            0.5 * spacing };

    // allocate grids for Tim's interpolator - the data is strided along dims
    std::unique_ptr< double[] > grid_ptr( new double[ n_dims * n_grid_points ] );

    // populate the grids
    for( int i = 0; i < n_dims; i++ )
    {
      for( int j = 0; j < n_grid_points; j++ )
      {
        grid_ptr[ i * n_grid_points + j ] = spacing * j + domain_start;
      }
    }

    double *grid_raw = grid_ptr.get();

    // calculate data lattice size
    int lattice_size = 1;

    for( int i = 0; i < n_dims; i++ )
    {
      lattice_size *= n_grid_points;
    }

    int lattice_cells = 1;
    for( int i = 0; i < n_dims; i++ )
    {
      lattice_cells *= ( n_grid_points - 1 );
    }

    std::vector< int > dims( n_dims, n_grid_points );

    // allocate lattice
    std::unique_ptr< double[] > lattice_ptr( new double[ lattice_size ] );

    double *lattice_data = lattice_ptr.get();

    // put lattice_data into a tensor object
    interpolated_field< double > 
    lattice( lattice_data, 
             dims.data(), 
             max_range_init.data(), 
             min_range_init.data(), 
             n_dims );

    // iterate lattice grid space with a round robin object
    round_robin_nd< n_dims, n_grid_points > lattice_indexer;

    round_robin_nd< n_dims, n_grid_points - 1 > cell_indexer;

    // create analytical function for testing
    auto analytical_function =
    []( std::array< double, n_dims > &real_space ) -> double
    {
      double d = 0;

      for( int j = 0; j < real_space.size(); j++ ) d += real_space[ j ] * real_space[ j ];

      d = 10 * std::cos( std::sqrt( d ) );

      return d;
    };

    /* performance checking loop, copy of above with other stuff pulled out  */
    auto const start = std::chrono::high_resolution_clock::now();
    for( int i = 0; i < lattice_cells; i++ )
    {
      cell_indexer.back_spin();

      auto lattice_gridpoint = cell_indexer.get_indices();

      std::array< double, n_dims > real_space;
      std::array< double, n_dims > real_space_offset;

      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      real_space[ j ] = lattice_gridpoint[ j ]*spacing+domain_start;

      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      {
        real_space_offset[ j ] = real_space[ j ] + offset[ j ];
      }

      double interpolated_value = lattice( real_space_offset.data() );
    }
    auto const stop = std::chrono::high_resolution_clock::now();

    double const duration_total =
    std::chrono::duration< double, std::chrono::nanoseconds::period >( stop - start ).count() / lattice_cells;

    /* performance checking loop, copy of above with other stuff pulled out  */
    auto const start_1 = std::chrono::high_resolution_clock::now();
    for( int i = 0; i < lattice_cells; i++ )
    {
      cell_indexer.back_spin();

      auto lattice_gridpoint = cell_indexer.get_indices();

      std::array< double, n_dims > real_space;
      std::array< double, n_dims > real_space_offset;

      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      real_space[ j ] = lattice_gridpoint[ j ]*spacing+domain_start;

      for( int j = 0; j < lattice_gridpoint.size(); j++ )
      {
        real_space_offset[ j ] = real_space[ j ] + offset[ j ];
      }

      double legacy_interpolation =
      interp3d(  
      real_space_offset[ 2 ],
      real_space_offset[ 1 ],
      real_space_offset[ 0 ],
      n_grid_points,
      n_grid_points,
      n_grid_points,
      grid_raw,
      grid_raw + n_grid_points,
      grid_raw + ( n_grid_points * 2 ),
      lattice_data );
    }

    auto const stop_1 = std::chrono::high_resolution_clock::now();

    double const duration_total_1 = 
    std::chrono::duration< double, std::chrono::nanoseconds::period >( stop_1 - start_1 ).count() / lattice_cells;

    std::cout << "legacy: " << duration_total_1 << " new: " << duration_total << std::endl;
  }
}















