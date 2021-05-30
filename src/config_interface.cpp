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
  lookup[ use::use_cuda ] = "USE_CUDA";
  lookup[ use::use_openmp ] = "USE_OPENMP";
  lookup[ use::use_mpi ] = "USE_MPI";
  lookup[ use::useionization ] = "USEIONIZATION";
  lookup[ use::use_ionization ] = "USE_IONIZATION";
  lookup[ use::userecombination ] = "USERECOMBINATION";
  lookup[ use::useperpdiffusion ] = "USEPERPDIFFUSION";
  lookup[ use::usecoulombcollisions ] = "USECOULOMBCOLLISIONS";
  lookup[ use::usefriction ] = "USEFRICTION";
  lookup[ use::useanglescattering ] = "USEANGLESCATTERING";
  lookup[ use::useheating ] = "USEHEATING";
  lookup[ use::usethermalforce ] = "USETHERMALFORCE";
  lookup[ use::usesurfacemodel ] = "USESURFACEMODEL";
  lookup[ use::usesheathefield ] = "USESHEATHEFIELD";
  lookup[ use::biased_surface ] = "BIASED_SURFACE";
  lookup[ use::usepresheathefield ] = "USEPRESHEATHEFIELD";
  lookup[ use::bfield_interp ] = "BFIELD_INTERP";
  lookup[ use::lc_interp ] = "LC_INTERP";
  lookup[ use::generate_lc ] = "GENERATE_LC";
  lookup[ use::efield_interp ] = "EFIELD_INTERP";
  lookup[ use::presheath_interp ] = "PRESHEATH_INTERP";
  lookup[ use::density_interp ] = "DENSITY_INTERP";
  lookup[ use::temp_interp ] = "TEMP_INTERP";
  lookup[ use::flowv_interp ] = "FLOWV_INTERP";
  lookup[ use::gradt_interp ] = "GRADT_INTERP";
  lookup[ use::odeint ] = "ODEINT";
  lookup[ use::fixedseeds ] = "FIXEDSEEDS";
  lookup[ use::fixed_seeds ] = "FIXED_SEEDS";
  lookup[ use::particleseeds  ] = "PARTICLESEEDS ";
  lookup[ use::geom_trace  ] = "GEOM_TRACE ";
  lookup[ use::geom_hash ] = "GEOM_HASH";
  lookup[ use::geom_hash_sheath ] = "GEOM_HASH_SHEATH";
  lookup[ use::particle_tracks ] = "PARTICLE_TRACKS";
  lookup[ use::particle_source_space ] = "PARTICLE_SOURCE_SPACE";
  lookup[ use::particle_source_energy ] = "PARTICLE_SOURCE_ENERGY";
  lookup[ use::particle_source_angle ] = "PARTICLE_SOURCE_ANGLE";
  lookup[ use::particle_source_file ] = "PARTICLE_SOURCE_FILE";
  lookup[ use::spectroscopy ] = "SPECTROSCOPY";
  lookup[ use::use3dtetgeom ] = "USE3DTETGEOM";
  lookup[ use::flux_ea ] = "FLUX_EA";
  lookup[ use::usecylsymm ] = "USECYLSYMM";
  lookup[ use::usefieldalignedvalues ] = "USEFIELDALIGNEDVALUES";
  lookup[ use::force_eval ] = "FORCE_EVAL";
  lookup[ use::check_compatibility ] = "CHECK_COMPATIBILITY";
  lookup[ use::use_sort ] = "USE_SORT";
  lookup[ use::use_adaptive_dt ] = "USE_ADAPTIVE_DT";
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
