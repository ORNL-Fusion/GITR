#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <iterator>

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
                      std::string module_path = "");

  /* get config sub-module */
  /* default behavior is to return the submodule itself */
  template< typename T = std::shared_ptr< class config_module_base > >
  T get( int key );

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
