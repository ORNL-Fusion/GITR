#pragma once
#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <iterator>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "netcdf.h"
#include "ncFile.h"
#include "ncVar.h"
#include "ncDim.h"
#include "libconfig.h++"

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

class netcdf_string_query
{
  public:

  netcdf_string_query( std::string netcdf_file_init = "" );

  ~netcdf_string_query();

  template< typename T >
  void operator()( std::string const &query_key, T &query_value ) const;

  template< typename T >
  void operator()( std::string const &query_key, T* query_value ) const;

  std::string netcdf_file_name = "empty";

  netCDF::NcFile nc;
};

template< typename T >
void netcdf_string_query::operator()( std::string const &query_key, T* query_value ) const
{
  std::cout << "Ahoy, ! Correct overload called for " << query_key << std::endl;

  netCDF::NcVar ncvar;

  ncvar = nc.getVar( query_key );

  if( ncvar.isNull() )
  {
    std::cout << "error: could not find " << query_key << std::endl;
  }

  else ncvar.getVar( query_value );
}

template< typename T >
void netcdf_string_query::operator()( std::string const &query_key, T &query_value ) const
{
  netCDF::NcDim nc_nx(nc.getDim(query_key));

  int n_x = nc_nx.getSize();

  query_value = n_x;
}


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
  void operator()( std::string const query_key, T &query_value ) const;

  template< typename T >
  void operator()( std::string const query_key, std::vector< T > &query_values ) const;

  private:

  libconfig::Config cfg;
};

template< typename T >
void libconfig_string_query::operator()( std::string const query_key,
                                         std::vector< T > &query_values ) const
{
  assert( query_values.size() == 0 );

  if( cfg.exists( query_key.c_str() ) == false )
  {
    throw 0;
  }

  auto const &setting = cfg.getRoot().lookup( query_key );

  /* this should be the local name of the setting */
  std::string setting_name( setting.getName() );

  /* check for array status */
  if( setting.isArray() )
  {
    int len = setting.getLength();

    for( int i = 0; i < len; i++ )
    {
      T value = setting[ i ];

      query_values.push_back( value );
    }
  }

  else
  {
    throw 0;
    //throw not_an_array( query_key );
  }
}

template< typename T >
void libconfig_string_query::operator()( std::string const query_key, T &query_value ) const
{
  if( cfg.exists( query_key.c_str() ) == false )
  {
    std::cout << "libconfig could not find " << query_key << ": " 
              << "leaving value at default: " << query_value
              << std::endl;
    return;
  }

  auto const &setting = cfg.getRoot().lookup( query_key );

  if( setting.isArray() )
  {
    //throw not_a_scalar( query_key );
    throw 0;
  }

  else
  {
    bool success = cfg.lookupValue( query_key, query_value );

    if( success == false ) 
    {
      //throw lookup_failed( query_key );
      std::cout << "libconfig " << query_key << " lookup failed!" << std::endl;
      throw 0;
    }
  }
}

class flags 
{
  public:

  int ionization = -1;
  int perp_diffusion = -1;
  int coulomb_collisions = -1;
  int thermal_force = -1;
  int surface_model = -1; 
  int sheath_efield = -1;
  int bfield_interp = -1;
  int efield_interp = -1;
  int density_interp = -1;
  int temp_interp = -1;
  int flowv_interp = -1;
  int gradt_interp = -1;
  int flux_ea = -1;
  int particle_source_file  = -1;
  int spectroscopy = -1;
  int use_3d_geom  = -1;
  int particle_tracks = -1;
  int fixed_seeds = -1;
  int geom_hash= -1;
  int geom_hash_sheath= -1;
  int cylsymm = -1;
  int force_eval = -1;
  int sort = -1;
  int adaptive_dt = -1;
  int surface_potential = -1;
  int particle_diagnostics = 0;
  int sheath_density = 0;
  int presheath_efield = 1;
  
  // Captain! Turn all this into a generate_flags() function and replace the old options
  CUDA_CALLABLE_MEMBER
  flags( libconfig_string_query const &query_metadata )
  {
    std::string const module_name = "flags";

    std::string sheath_density_str = "USE_SHEATH_DENSITY";
    query_metadata( module_name + "." + sheath_density_str, sheath_density );

    std::string particle_diagnostics_str = "USE_PARTICLE_DIAGNOSTICS";
    query_metadata( module_name + "." + particle_diagnostics_str, particle_diagnostics );

    std::string ionization_str = "USE_IONIZATION";
    query_metadata( module_name + "." + ionization_str, ionization );

    std::string perp_diffusion_str = "USEPERPDIFFUSION";
    query_metadata( module_name + "." + perp_diffusion_str, perp_diffusion );

    std::string coulomb_collisions_str = "USECOULOMBCOLLISIONS";
    query_metadata( module_name + "." + coulomb_collisions_str, coulomb_collisions );

    std::string thermal_force_str = "USETHERMALFORCE";
    query_metadata( module_name + "." + thermal_force_str, thermal_force );

    std::string surface_model_str = "USESURFACEMODEL";
    query_metadata( module_name + "." + surface_model_str, surface_model );

    std::string sheath_efield_str = "USESHEATHEFIELD";
    query_metadata( module_name + "." + sheath_efield_str, sheath_efield );

    std::string bfield_interp_str = "BFIELD_INTERP";
    query_metadata( module_name + "." + bfield_interp_str, bfield_interp );

    std::string efield_interp_str =  "EFIELD_INTERP";
    query_metadata( module_name + "." + efield_interp_str,  efield_interp );

    std::string density_interp_str = "DENSITY_INTERP";
    query_metadata( module_name + "." + density_interp_str ,  density_interp );

    std::string temp_interp_str = "TEMP_INTERP";
    query_metadata( module_name + "." + temp_interp_str ,  temp_interp );

    std::string flowv_interp_str = "FLOWV_INTERP";
    query_metadata( module_name + "." + flowv_interp_str ,  flowv_interp );

    std::string gradt_interp_str = "GRADT_INTERP";
    query_metadata( module_name + "." + gradt_interp_str ,  gradt_interp );

    std::string particle_tracks_str = "PARTICLE_TRACKS";
    query_metadata( module_name + "." + particle_tracks_str , particle_tracks );

    std::string particle_source_file_str  = "PARTICLE_SOURCE_FILE";
    query_metadata( module_name + "." + particle_source_file_str  , particle_source_file  );

    std::string spectroscopy_str  = "SPECTROSCOPY";
    query_metadata( module_name + "." + spectroscopy_str  , spectroscopy );

    std::string use_3d_geom_str  = "USE3DTETGEOM";
    query_metadata( module_name + "." + use_3d_geom_str  , use_3d_geom  );

    std::string flux_ea_str = "FLUX_EA";
    query_metadata( module_name + "." + flux_ea_str , flux_ea  );

    std::string fixed_seeds_str = "FIXED_SEEDS";
    query_metadata( module_name + "." + fixed_seeds_str , fixed_seeds );

    std::string geom_hash_str = "GEOM_HASH";
    query_metadata( module_name + "." + geom_hash_str , geom_hash );

    std::string geom_hash_sheath_str = "GEOM_HASH_SHEATH";
    query_metadata( module_name + "." + geom_hash_sheath_str , geom_hash_sheath );

    std::string cylsymm_str = "USECYLSYMM";
    query_metadata( module_name + "." + cylsymm_str , cylsymm );
    
    std::string force_eval_str  = "FORCE_EVAL";
    query_metadata( module_name + "." + force_eval_str  , force_eval  );

    std::string sort_str  = "USE_SORT";
    query_metadata( module_name + "." + sort_str  , sort  );

    std::string adaptive_dt_str = "USE_ADAPTIVE_DT";
    query_metadata( module_name + "." + adaptive_dt_str , adaptive_dt );

    std::string surface_potential_str = "USE_SURFACE_POTENTIAL";
    query_metadata( module_name + "." + surface_potential_str, surface_potential );
  }
};
