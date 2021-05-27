#include "config_interface.h"

libconfig_string_query::libconfig_string_query( std::string libconfig_file )
{
  /* open the file */
  try
  {
    cfg.readFile("example.cfg");
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

impurity_particle_source::
impurity_particle_source( std::string const module_path )
  :
  config_module_base( module_path )
{
}

template< typename T >
void impurity_particle_source::emulate_template( int key, T &val )
{
  auto access = lookup.find( key );
  if( access == lookup.end() )
  {
    std::cout << "error: value key not found" << std::endl;
    exit(0);
  }

  query( module_path() + access.second, val );
}

/* explicit instantiations of that template */
template void impurity_particle_source::emulate_template<int>( int key, int &val );
template void impurity_particle_source::emulate_template<float>( int key, float &val );
template void impurity_particle_source::emulate_template<double>( int key, double &val );
template void 
impurity_particle_source::emulate_template<std::string>( int key, std::string &val );


/* get a config submodule */
std::shared_ptr< config_module_base >
impurity_particle_source::get( int key )
{
  auto access = sub_modules.find( key );

  if( access == lookup.end() )
  {
    std::cout << "error: config_module key not found" << std::endl;
    exit( 0 );
  }

  return access.second;
}

/* These are not templated because they are overrides  */
void impurity_particle_source::get( int key, int &val )
{
  emulate_template( key, val );
}

void impurity_particle_source::get( int key, float &val )
{
  emulate_template( key, val );
}

void impurity_particle_source::get( int key, double &val )
{
  emulate_template( key, val );
}

void impurity_particle_source::get( int key, std::string &val )
{
  emulate_template( key, val );
}
