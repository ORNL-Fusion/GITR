#include <iostream>
#include "test_utils.hpp"
#include "config_interface.h"

/* auto-generated header contains filepath of unit testing files
   defines LIBCONFIG_UNIT_TEST_FILE */
#include "test_data_filepath.hpp"

TEST_CASE( "Simulation Configuration - not fully implemented" )
{
  /* Open the file and get a field out */
  class libconfig_string_query query( CONFIG_INTERFACE_UNIT_TEST_FILE );

  SECTION( "read raw libconfig data" )
  {
    std::string search = "surfaces.flux.E";

    double e = 0;
    
    query( search, e );

    REQUIRE( e == 1000 );
  }

  /* build a config module similar to how you would in a real situation and test it */
  /* Just get this to compile and get the skeleton set up. Commit that work and
     then you can stop for today */
  SECTION( "config_module" )
  {
    /* constructor should create all the submodules necessary underneath itself */
    auto impurity_particle_source = std::make_shared< class impurity_particle_source >( query );

    /* get a value from this config module */
    int source_material_z = 
    impurity_particle_source->get<int>( impurity_particle_source::source_material_z );

    REQUIRE( source_material_z == 13 );

    /* get a child module and obtain a value from it */
    auto ionization = impurity_particle_source->get( impurity_particle_source::ionization );

    auto recombination =
    impurity_particle_source->get( impurity_particle_source::recombination );

    std::string dense_grid_string =
    ionization->get<std::string>( ionization_process::dense_grid_string );

    REQUIRE( dense_grid_string == "n_Densities_Ionize" );

    REQUIRE( dense_grid_string == impurity_particle_source
             ->get( impurity_particle_source::ionization )
             ->get< std::string >( ionization_process::dense_grid_string ) );

    dense_grid_string =
    recombination->get<std::string>( ionization_process::dense_grid_string );

    REQUIRE( dense_grid_string == "n_Densities_Recombine" );

    REQUIRE( dense_grid_string == impurity_particle_source
             ->get( impurity_particle_source::recombination )
             ->get< std::string >( ionization_process::dense_grid_string ) );

    class use use( query );

    int spectroscopy = use.get<int>( use::spectroscopy );

    REQUIRE( spectroscopy == 3);
  }

  /* Test the exception handling - also create a malformed class */
  SECTION( "misc" )
  {
    auto geometry = std::make_shared< class geometry >( query );

    /* get an array value from this config module */
    auto slope = geometry->get< std::vector< double > >( geometry::slope );

    REQUIRE( slope == std::vector< double >{ 1e12, 1e12 } );

    /* try opening a bad file too. After this, be done with GITR this week please */

    /* Captain! Catch the exceptions it throws - expect a specific key */

    /* also handle the case where the desired config option is not available */

    /* the key was not registered with the lookup case, the key does not exist in the 
       libconfig file */

    /* What else needs to be handled? */
    /*

       Config options do need to be cached... but if you are accessing them via
       ->get
       ->get
       ...
       successive pointer dereferences, you're going to have a lot of issues with performance.
       Is it possible to avoid dereferencing long chains of them in the functors?

       Just leave a comment and don't do that non-performant behavior
    
    */
  }
}












