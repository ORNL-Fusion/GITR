#include "config_interface.h"

libconfig_string_query::libconfig_string_query( std::string libconfig_file )
{
  /* open the file */
  try
  {
    cfg.readFile( libconfig_file.c_str() );
  }

  /* bad file */
  catch(const libconfig::FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    exit( 0 );
  }

  /* bad format */
  catch(const libconfig::ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    exit( 0 );
  }
}

geometry::
geometry( class libconfig_string_query const &query,
          std::string module_path )
  :
  config_module_base( query, module_path )
{
  lookup[ geometry::slope ] = "slope";
  lookup[ geometry::intercept ] = "intercept";
  lookup[ geometry::length ] = "length";
  lookup[ geometry::z ] = "Z";
  lookup[ geometry::surface ] = "surface";
  lookup[ geometry::in_dir ] = "in_dir";
  lookup[ geometry::periodic ] = "periodic";
  lookup[ geometry::x1 ] = "x1";
  lookup[ geometry::x2 ] = "x2";
  lookup[ geometry::y1 ] = "y1";
  lookup[ geometry::y2 ] = "y2";
  lookup[ geometry::z1 ] = "z1";
  lookup[ geometry::z2 ] = "z2";
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
  lookup[ impurity_particle_source::recombination ] = "recombination";

  /* create ionization config submodule */
  generate_sub_module< ionization_process >
  ( impurity_particle_source::ionization );

  /* create recombination config submodule */
  generate_sub_module< ionization_process >
  ( impurity_particle_source::recombination );
}

template< typename T >
void config_module_base::generate_sub_module( int key )
{
  std::shared_ptr< T >
  sub_module(
  new T( query, 
         get_module_path() +
         "." +
         lookup[ key ] ) );

  sub_modules[ key ] = sub_module;
}

ionization_process::
  ionization_process( class libconfig_string_query const &query,
                      std::string module_path )
  :
  config_module_base( query, module_path )
{ 
  lookup[ ionization_process::file_string ] = "fileString";
  lookup[ ionization_process::temp_grid_string ] = "TempGridString";
  lookup[ ionization_process::dense_grid_string ] = "DensGridString";
  lookup[ ionization_process::charge_state_string ] = "nChargeStateString";
  lookup[ ionization_process::temp_grid_var_name ] = "TempGridVarName";
  lookup[ ionization_process::dense_grid_var_name ] = "DensGridVarName";
  lookup[ ionization_process::coeff_var_name ] = "CoeffVarName";
}

use::use( class libconfig_string_query const &query,
          std::string module_path )
  :
  config_module_base( query, module_path )
{
  lookup[ use::cuda ] = "USE_CUDA";
  lookup[ use::use_openmp ] = "USE_OPENMP";
  lookup[ use::mpi ] = "USE_MPI";
  lookup[ use::ionization ] = "USE_IONIZATION";
  lookup[ use::perp_diffusion ] = "USEPERPDIFFUSION";
  lookup[ use::coulomb_collisions ] = "USECOULOMBCOLLISIONS";
  //lookup[ use::friction ] = "USEFRICTION";
  //lookup[ use::angle_scattering ] = "USEANGLESCATTERING";
  //lookup[ use::heating ] = "USEHEATING";
  lookup[ use::thermal_force ] = "USETHERMALFORCE";
  lookup[ use::surface_model ] = "USESURFACEMODEL";
  lookup[ use::sheath_efield ] = "USESHEATHEFIELD";
  // hardcoded to 0 for now
  //lookup[ use::biased_surface ] = "BIASED_SURFACE";
  //lookup[ use::presheath_efield ] = "USEPRESHEATHEFIELD";
  lookup[ use::bfield_interp ] = "BFIELD_INTERP";
  lookup[ use::efield_interp ] = "EFIELD_INTERP";
  lookup[ use::presheath_interp ] = "PRESHEATH_INTERP";
  lookup[ use::density_interp ] = "DENSITY_INTERP";
  lookup[ use::temp_interp ] = "TEMP_INTERP";
  lookup[ use::flowv_interp ] = "FLOWV_INTERP";
  lookup[ use::gradt_interp ] = "GRADT_INTERP";
  lookup[ use::odeint ] = "ODEINT";
  lookup[ use::fixed_seeds ] = "FIXED_SEEDS";
  lookup[ use::geom_trace  ] = "GEOM_TRACE ";
  lookup[ use::geom_hash ] = "GEOM_HASH";
  lookup[ use::geom_hash_sheath ] = "GEOM_HASH_SHEATH";
  lookup[ use::particle_tracks ] = "PARTICLE_TRACKS";
  lookup[ use::particle_source_space ] = "PARTICLE_SOURCE_SPACE";
  lookup[ use::particle_source_energy ] = "PARTICLE_SOURCE_ENERGY";
  lookup[ use::particle_source_angle ] = "PARTICLE_SOURCE_ANGLE";
  lookup[ use::particle_source_file ] = "PARTICLE_SOURCE_FILE";
  lookup[ use::spectroscopy ] = "SPECTROSCOPY";
  lookup[ use::use_3d_geom ] = "USE3DTETGEOM";
  lookup[ use::flux_ea ] = "FLUX_EA";
  lookup[ use::cylsymm ] = "USECYLSYMM";
  lookup[ use::field_aligned_values ] = "USEFIELDALIGNEDVALUES";
  lookup[ use::force_eval ] = "FORCE_EVAL";
  //lookup[ use::compatibility_check ] = "CHECK_COMPATIBILITY";
  lookup[ use::sort ] = "USE_SORT";
  lookup[ use::adaptive_dt ] = "USE_ADAPTIVE_DT";
  lookup[ use::surface_potential ] = "USE_SURFACE_POTENTIAL";
  lookup[ use::sheath_model_type ] = "SHEATH_MODEL_TYPE";
    lookup[ use::nspecies ] = "NSPECIES";
}

config_module_base::config_module_base( class libconfig_string_query const &query,
                                        std::string module_path )
  :
  module_path( module_path ),
  query( query )
{ }

/* general "get" for values  */
template< typename T >
T config_module_base::get( int key )
{
  auto access = lookup.find( key );

  if( access == lookup.end() )
  {
    throw unregistered_config_mapping( key );
  }

  T val;

  try
  {
    query( get_module_path() + "." + access->second, val );
  }

  catch( class invalid_key const &exception )
  {
    std::string const error_message{ exception.what() };

    std::string const error_key{ exception.get_key() };

    std::cout << error_message << error_key << std::endl;

    throw invalid_key( error_key );
  }

  catch( class lookup_failed const &exception )
  {
    std::string const error_message{ exception.what() };

    std::string const error_key{ exception.get_key() };

    std::cout << error_message << error_key << std::endl;

    throw lookup_failed( error_key );
  }
  
  return val;
}

/* instantiations for singleton config setting values */
template int config_module_base::get<int>( int key );
template float config_module_base::get<float>( int key );
template double config_module_base::get<double>( int key );
template bool config_module_base::get<bool>( int key );
template std::string config_module_base::get<std::string>( int key ); 

/* instantiations for vector/array config setting values */
template std::vector< int > config_module_base::get< std::vector< int > >( int key );
template std::vector< float >
config_module_base::get< std::vector< float > >( int key );
template std::vector< double >
config_module_base::get< std::vector< double > >( int key );
template std::vector< bool > config_module_base::get< std::vector< bool > >( int key );
template std::vector< std::string >
config_module_base::get< std::vector< std::string > >( int key ); 

/* template specialization for non-specified type: specialization for default T from header */
template<>
std::shared_ptr< config_module_base >
config_module_base::get( int key )
{
  auto access = sub_modules.find( key );

  if( access == sub_modules.end() )
  {
    throw( unregistered_config_mapping( key ) );
  }

  return access->second;
}
