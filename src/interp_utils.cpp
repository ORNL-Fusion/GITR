#include "interpolator.h"
#include "interp_utils.h"
#include "operation.hpp"

domain::domain( double domain_start_in,
        double domain_end_in,
        int const n_grid_points_in,
        int const n_dims_in )
  :
  domain_start( domain_start_in ),
  domain_end( domain_end_in ),
  n_grid_points( n_grid_points_in ),
  n_dims( n_dims_in ),
  dims( new int[ n_dims ] ),
  min_range_init( new double[ n_dims ] ),
  max_range_init( new double[ n_dims ] ),
  spacing( ( domain_end - domain_start ) / ( n_grid_points - 1 ) )
{ 
  for( int i = 0; i < n_dims; i++ ) dims[ i ] = n_grid_points;

  for( int i = 0; i < n_dims; i++ ) min_range_init[ i ] = domain_start;

  for( int i = 0; i < n_dims; i++ ) max_range_init[ i ] = domain_end;

  // offset + cell_corner_point = cell_midpoint
  offset = { 0.5 * spacing, 0.5 * spacing, 0.5 * spacing };

  lattice_size = 1;

  for( int i = 0; i < n_dims; i++ )
  {
    lattice_size *= n_grid_points;
  }

  lattice_cells = 1;

  for( int i = 0; i < n_dims; i++ )
  {
    lattice_cells *= ( n_grid_points - 1 );
  }
}

legacy_domain::legacy_domain( domain &d )
  :
  grid_ptr( new double[ d.n_dims * d.n_grid_points ] ),
  grid_raw( grid_ptr.get() )
{
  // populate square grids
  // data is strided along dims
  for( int i = 0; i < d.n_dims; i++ )
  {
    for( int j = 0; j < d.n_grid_points; j++ )
    {
      grid_ptr[ i * d.n_grid_points + j ] = d.spacing * j + d.domain_start;
    }
  }
}

std::array< double, 3 > domain::generate_point( int index )
{
  int r = index;
  std::array< int, 3 > lattice_gridpoint{ 0, 0, 0 };
  for( int j = 0; j < n_dims; j++ )
  {
    int m = 1;
    for( int k = 0; k < n_dims - j - 1; k++ )
    {
      m *= ( n_grid_points - 1 );
    }

    int a = r / m;

    r -= ( a * m );

    lattice_gridpoint[ j ] = a;
  }

  /* convert cell coordinate into xyz realspace coordinate for bottom corner of cell */
  std::array< double, 3 > real_space;
  for( int j = 0; j < lattice_gridpoint.size(); j++ )
  real_space[ j ] = lattice_gridpoint[ j ]*spacing+domain_start;

  /* convert realspace coordinate into the center of the cell instead of the corner */
  std::array< double, 3 > real_space_offset;
  for( int j = 0; j < lattice_gridpoint.size(); j++ )
  {
    real_space_offset[ j ] = real_space[ j ] + offset[ j ];
  }

  return real_space_offset;
}

populated_lattice::populated_lattice( domain const &d, 
                                      std::function<double(std::array< double, 3 >&)> 
                                      analytical_function ) 
  :
  lattice_ptr( new double[ d.lattice_size ] ),
  ifield( lattice_ptr.get(), 
      d.dims.get(), 
      d.max_range_init.get(), 
      d.min_range_init.get(), 
      d.n_dims )
{
  // populate ifield with data - Captain! Make sure this works then take out the indexer
  //round_robin_nd< d.n_dims, d.n_grid_points > lattice_indexer;

  for( int i = 0; i < d.lattice_size; i++ )
  {
    /* grid-space */
    //lattice_indexer.back_spin();
    //auto lattice_gridpoint_0 = lattice_indexer.get_indices();

    // iterate over the coordinates
    // z y x order
    int r = i;
    std::array< int, 3 > lattice_gridpoint{ 0, 0, 0 };
    for( int j = 0; j < d.n_dims; j++ )
    {
      int m = 1;
      for( int k = 0; k < d.n_dims - j - 1; k++ )
      {
        m *= d.n_grid_points;
      }

      int a = r / m;

      r -= ( a * m );

      lattice_gridpoint[ j ] = a;
    }

    //if( lattice_gridpoint_0 != lattice_gridpoint ) exit( 0 );
    std::array< double, 3 > real_space;

    for( int j = 0; j < lattice_gridpoint.size(); j++ )
      real_space[ j ] = lattice_gridpoint[ j ]*d.spacing+d.domain_start;

    double analytical_value = analytical_function( real_space );

    ifield.set( analytical_value, lattice_gridpoint.data() );

    double stored_analytical_value = ifield.get( lattice_gridpoint.data() );

    if( analytical_value != stored_analytical_value ) exit( 0 );
    }
}

output_builder::output_builder( std::filesystem::path &output_file_name, domain &d_in )
  :
  data_stack_3d( output_file_name ),
  d( d_in )
{
}

double output_builder::get_mean_difference()
{
  return total_difference / d.lattice_cells;
}

double output_builder::get_mean_legacy_difference()
{
  return total_legacy_difference / d.lattice_cells;
}

void output_builder::process_timestep( std::array< double, 3 > &real_space_offset,
                       double analytical_offset_value, 
                       double legacy_interpolation,
                       double interpolated_value,
                       int timestep )
{
  // stack the timestep info into output_file_name, verify it is identical
  // print output to screen
  // calculate total and legacy difference
  // sweep the gridpoints to get multiple values for this...

  /* output results */
  /* take a projection */
  if( real_space_offset[ 1 ] == d.offset[ 1 ] )
  {
    // columns: cell_index_i x y z f_analytical f_old f_new
    std::array< double, 6 > 
      d{ real_space_offset[ 2 ], real_space_offset[ 1 ], real_space_offset[ 0 ],
        analytical_offset_value, legacy_interpolation, interpolated_value };

    std::span< double, 6 > data_span( d );

    data_stack_3d.stack( timestep, data_span );
  }

  double difference = 
    std::abs( ( interpolated_value - analytical_offset_value ) / analytical_offset_value );

  total_difference += difference;

  double legacy_difference =
    std::abs( ( legacy_interpolation - analytical_offset_value ) / analytical_offset_value );

  total_legacy_difference += legacy_difference;

  /* print test code - check them */
  if( timestep % 1000000  == 0 )
  {
    std::cout << timestep << "/" << d.lattice_cells << " = " 
      << (double)timestep / (double)d.lattice_cells * 100 << std::endl;

    std::cout << "new method difference: " << std::setprecision( 10 )
      << difference
      << std::endl;

    std::cout << "legacy method difference: " << std::setprecision( 10 )
      << legacy_difference
      << std::endl;

    std::cout << std::setprecision( 10 ) << "discrepancy: "
      << ( interpolated_value - legacy_interpolation ) / analytical_offset_value
      << std::endl;
  }
}



























