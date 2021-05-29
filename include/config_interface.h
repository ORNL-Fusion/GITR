#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <iterator>

#include "libconfig.h++"

/* instead of yagni, make an extendability note */

/* This class wraps the libconfig object and inserts
   a key-string mapping in between. This creates an
   internal separation between the internal representation
   of the configuration and the I/O details of obtaining it
   */

/* extensibility note: augment this to derive from a string query base class */
class libconfig_string_query
{
  public:

  libconfig_string_query( std::string libconfig_file = "" );

  template< typename T >
  void operator()( std::string const query_key, T &query_value ) const
  {
    bool success = cfg.lookupValue( query_key, query_value );

    if( success == false ) 
    {
      std::cout << "invalid query key: " << query_key << std::endl;
      exit(0);
    }
  }

  private:

  libconfig::Config cfg;
};

class config_module_base
{
  public:

  config_module_base( class libconfig_string_query const &query,
                      std::string module_path = "" );

  /* get config sub-module */
  std::shared_ptr< config_module_base > get( int key );

  /* a config value */
  template< typename T >
  void get( int key, T &val );

  protected:

  std::string const &get_module_path() { return module_path; }

  /* data structure maps internal representation ---> string name of config in config file */
  std::unordered_map< int, std::string > lookup;

  /* data structure maps internal representation ---> string name of a child module in 
     config file */
  std::unordered_map< int, std::shared_ptr< config_module_base > > sub_modules;

  /* string name of this modules path in the config file */
  std::string const module_path;

  /* I/O interface to the config file */
  class libconfig_string_query const &query;
};

/* Captain! Move the derived stuff into different files! Then you can implement these
   specifically for the unit tests and migrate them over to the main codebase if
   they end up being useful. Also add documentation about how to add new configurations */
class ionization_process final : public config_module_base
{
  public:

  /* Captain! Document this */
  enum : int
  {
    file_string,
    temp_grid_string,
    dense_grid_string,
    charge_state_string,
    temp_grid_var_name,
    dense_grid_var_name,
    coeff_var_name
  };

  ionization_process( class libconfig_string_query const &query,
                      std::string module_path );

  ionization_process( std::string module_path = 
  "this module represents several instances" );
};

class impurity_particle_source final : public config_module_base
{
  public:

  /* Captain! Document this */
  enum : int
  {
    source_material_z,
    ionization
  };

  impurity_particle_source( class libconfig_string_query const &query,
                      std::string module_path = "impurityParticleSource" );
};

/* top-level config module */
class gitr_config_interface final : public config_module_base
{
  public:

  gitr_config_interface();
};
