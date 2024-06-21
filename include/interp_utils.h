#include <functional>
#include "io.h"

class domain
{
  public:

  domain( double domain_start_in,
          double domain_end_in,
          int const n_grid_points_in,
          int const n_dims_in );

  std::array< double, 3 > generate_point( int index );

  double const domain_start;
  double const domain_end;
  int const n_grid_points;
  int const n_dims;
  std::unique_ptr< int[] > dims;
  std::unique_ptr< double[] > min_range_init;
  std::unique_ptr< double[] > max_range_init;
  double const spacing;
  std::array< double, 3 > offset;
  int lattice_size;
  int lattice_cells;
};

class legacy_domain
{
  public:

  legacy_domain( domain &d );

  std::unique_ptr< double[] > grid_ptr;

  double *grid_raw;
};

class populated_lattice
{
  public:

  populated_lattice( domain const &d, 
                     std::function<double(std::array< double, 3 >&)> 
                     analytical_function ); 

  std::unique_ptr< double[] > lattice_ptr;

  interpolated_field< double > ifield;
};

class output_builder
{
  public:

  output_builder( std::filesystem::path &output_file_name, domain &d_in );

  void process_timestep( std::array< double, 3 > &real_space_offset,
                         double analytical_offset_value, 
                         double legacy_interpolation,
                         double interpolated_value,
                         int timestep );

  double get_mean_difference();
  double get_mean_legacy_difference();

  csv_row_stacker< 6 > data_stack_3d;
  domain &d;

  double total_difference{ 0 };
  double total_legacy_difference{ 0 };
};















