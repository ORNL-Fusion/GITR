#pragma once
#include "config_interface.h"
#include "array.h"
#include <cstdlib>

class temperature_grid_type
{
  public:

  temperature_grid_type( libconfig_string_query const &query_metadata, 
                         std::string parent_module_name )
    :
    full_module_path( parent_module_name + "." + module_name + "." )
  {
    // metadata keys (put in a namespace later?)
    std::string raw_input_file_key = "fileString";

    // meta keys have values that serve as keys to a final value
    std::string r_grid_size_meta_key = "gridNrString";
    std::string y_grid_size_meta_key = "gridNyString";
    std::string z_grid_size_meta_key = "gridNzString";

    std::string r_grid_meta_key = "gridRString";
    std::string y_grid_meta_key = "gridYString";
    std::string z_grid_meta_key = "gridZString";

    std::string ion_temp_data_meta_key = "IonTempString";
    std::string electron_temp_data_meta_key = "ElectronTempString";
    std::string ion_temp_data_key;
    std::string electron_temp_data_key;

    // create netcdf file reader object:
    std::string raw_input_file_path;

    query_metadata( full_module_path + raw_input_file_key, raw_input_file_path );

    // step 0: meta_key --> key
    std::string r_grid_size_key;
    std::string z_grid_size_key;
    query_metadata( full_module_path + r_grid_size_meta_key, r_grid_size_key );
    query_metadata( full_module_path + z_grid_size_meta_key, z_grid_size_key );

    query_metadata( full_module_path + ion_temp_data_meta_key, ion_temp_data_key );
    query_metadata( full_module_path + electron_temp_data_meta_key, electron_temp_data_key );

    // query query_raw_input using the keys and values above
    netcdf_string_query const query_raw_input( "input/" + raw_input_file_path );

    int r_grid_size = -1;
    int z_grid_size = -1;
    query_raw_input( r_grid_size_key, r_grid_size );
    query_raw_input( z_grid_size_key, z_grid_size );

    // check grid sizes here:
    if( r_grid_size == -1 ) std::cout << "error: r_grid_size is -1";
    if( z_grid_size == -1 ) std::cout << "error: z_grid_size is -1";

    ti.resize( r_grid_size * z_grid_size );
    te.resize( r_grid_size * z_grid_size );

    query_raw_input( ion_temp_data_key, ti.data() );
    query_raw_input( electron_temp_data_key, te.data() );
  }

  // members:
  sim::Array< double > ti;
  sim::Array< double > te;

  // include these as sim::array initially and then pass to tensor class
  // ti data
  // te data

  /*
     Temperature =
     {
      #ti = 10.0;
      #te = 10.0;
      fileString = "profiles.nc";
      gridNrString = "nR";
      gridNzString = "nZ";
      gridRString = "gridR";
      gridZString = "gridZ";
      IonTempString = "ti";
      ElectronTempString = "te";
      }
   */
  std::string const module_name = "Temperature";
  std::string const full_module_path;
};

class background_plasma
{
  public:

  std::string const module_name = "backgroundPlasmaProfiles";
  
  temperature_grid_type temperature_grid;

  background_plasma( libconfig_string_query const &query_metadata )
    :
    temperature_grid( query_metadata, module_name )
  { }
};
