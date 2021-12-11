#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <iterator>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "config_interface_exceptions.h"
#include "libconfig.h++"

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
  void operator()( std::string const &query_key, T &query_value ) const;

  template< typename T >
  void operator()( std::string const &query_key, std::vector< T > &query_values ) const;

  private:

  libconfig::Config cfg;
};

template< typename T >
void libconfig_string_query::operator()( std::string const &query_key,
                                         std::vector< T > &query_values ) const
{
  assert( query_values.size() == 0 );

  if( cfg.exists( query_key.c_str() ) == false )
  {
    throw invalid_key(query_key);
  }

  auto const &setting = cfg.getRoot().lookup( query_key );

  /* this should be the local name of the setting */
  std::string setting_name( setting.getName() );

  std::cout << "Captain! setting name: " << setting_name << std::endl;
  std::cout << "Captain! query key: " << query_key << std::endl;

  /* Ahoy, Captain! I think the issue here is that the "query_key" is geom.slope and the setting
     name is "slope". setting_name should be "geom" and query_key should be "slope" */
  /* How can we fix this? I think there is a way to get the path not the name */
  /* Ahoy, Captain! Use getParent and getPath and all that */

  /* check for array status */
  if( setting.isArray() )
  {
    int len = setting.getLength();

    for( int i = 0; i < len; i++ )
    {
      T value = setting[ i ];
      std::cout << "Captain! setting: " << value << std::endl;
      if( setting[ i ].lookupValue( setting_name.c_str(), value ) == false )
      {
        std::cout << "Captain! Type is: " << typeid( T ).name() << std::endl;
      }
      query_values.push_back( value );
    }
  }

  else
  {
    std::cout << "Ahoy! This setting is not a vector. Create a new exception for it"
              << std::endl;

    exit( 0 );
  }
}

template< typename T >
void libconfig_string_query::operator()( std::string const &query_key, T &query_value ) const
{
  if( cfg.exists( query_key.c_str() ) == false )
  {
    throw( invalid_key( query_key ) );
  }

  else
  {
    bool success = cfg.lookupValue( query_key, query_value );

    if( success == false ) 
    {
      std::cout << "Ahoy! lookup failed on: " << query_key << std::endl;
      throw( lookup_failed( query_key ) );
    }
  }
}

class config_module_base
{
  public:

  config_module_base( class libconfig_string_query const &query,
                      std::string module_path = "");

  /* I literally want "T" itself to be a template template parameter:
     I want T to mean std::vector< T >*/
  template< typename T = std::shared_ptr< class config_module_base > >
  T get( int key );

  /* Captain! To fix this, you need to create a "get" that calls the query operator on
     vectors instead of scalars. Rename this to get_vector() */
  /*
  template< typename T >
  std::vector< T > get( int key );
  */

  template< typename T >
  void generate_sub_module( int key );

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

class geometry final : public config_module_base
{
  public:

  enum : int
  {
    slope,
    intercept,
    length,
    z,
    surface,
    in_dir,
    periodic,
    x1,
    x2,
    y1,
    y2,
    z1,
    z2
  };

  geometry( class libconfig_string_query const &query,
            std::string module_path = "geom" );
};

class ionization_process final : public config_module_base
{
  public:

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
};

class impurity_particle_source final : public config_module_base
{
  public:

  enum : int
  {
    source_material_z, //this is the material in z
    ionization,
    recombination
  };

  impurity_particle_source( class libconfig_string_query const &query,
                            std::string module_path = "impurityParticleSource" );
};

class use final : public config_module_base
{
  public:

  enum : int
  {
    cuda,
    use_openmp,
    mpi,
    useionization,
    use_ionization,
    userecombination,
    useperpdiffusion,
    usecoulombcollisions,
    usefriction,
    useanglescattering,
    useheating,
    usethermalforce,
    usesurfacemodel,
    usesheathefield,
    biased_surface,
    usepresheathefield,
    bfield_interp,
    lc_interp,
    generate_lc,
    efield_interp,
    presheath_interp,
    density_interp,
    temp_interp,
    flowv_interp,
    gradt_interp,
    odeint,
    fixedseeds,
    fixed_seeds,
    particleseeds ,
    geom_trace ,
    geom_hash,
    geom_hash_sheath,
    particle_tracks,
    particle_source_space,
    particle_source_energy,
    particle_source_angle,
    particle_source_file,
    spectroscopy,
    use3dtetgeom,
    flux_ea,
    usecylsymm,
    usefieldalignedvalues,
    force_eval,
    compatibility_check,
    use_sort,
    use_adaptive_dt
  };

  use( class libconfig_string_query const &query,
       std::string module_path = "flags" );
};

/* top-level config module */
class gitr_config_interface final : public config_module_base
{
  public:

  gitr_config_interface();
};
