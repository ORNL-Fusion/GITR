#include "config_interface.h"

libconfig_string_query::libconfig_string_query( std::string libconfig_file )
{
  /* open the file */
  try
  {
    cfg.readFile( libconfig_file );
  }

  catch(const libconfig::FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
  }

  catch(const libconfig::ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
  }
}

/* get a config submodule */
std::shared_ptr< config_module_base >
config_module_base::get( int key )
{
  auto access = sub_modules.find( key );

  if( access == sub_modules.end() )
  {
    std::cout << "error: config_module key not found" << std::endl;
    exit( 0 );
  }

  return access->second;
}

/* geta config value */
template< typename T >
void config_module_base::get( int key, T &val )
{
  auto access = lookup.find( key );
  if( access == lookup.end() )
  {
    std::cout << "error: value key not found" << std::endl;
    exit(0);
  }

  query( get_module_path() + "." + access->second, val );
}

impurity_particle_source::
impurity_particle_source( class libconfig_string_query const &query,
                          std::string module_path )
  :
  config_module_base( query, module_path )
{ 
  /* declare values */
  lookup[ impurity_particle_source::source_material_z ] = "source_material_Z";
  lookup[ impurity_particle_source::ionization ] = "ionization";

  /* create ionization submodule */
  std::shared_ptr< ionization_process >
  ionization(
  new ionization_process( query, 
                          get_module_path() +
                          "." +
                          lookup[ impurity_particle_source::ionization ] ) );

  sub_modules[ impurity_particle_source::ionization ] = ionization;
}

ionization_process::
  ionization_process( class libconfig_string_query const &query,
                      std::string module_path )
  :
  config_module_base( query, module_path )
{ 
  lookup[ ionization_process::dense_grid_string ] = "DensGridString";
}

config_module_base::config_module_base( class libconfig_string_query const &query,
                                        std::string module_path )
  :
  module_path( module_path ),
  query( query )
{ }

/* explicit instantiations of that template */
template void config_module_base::get<int>( int key, int &val );
template void config_module_base::get<float>( int key, float &val );
template void config_module_base::get<double>( int key, double &val );
template void config_module_base::get<std::string>( int key, std::string &val ); 
