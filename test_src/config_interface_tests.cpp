#include <iostream>
#include "test_utils.hpp"
#include "config_interface.h"
/* contains filepath of unit testing files */
#include "test_data_filepath.hpp"

TEST_CASE( "Simulation Configuration" )
{
  /* Open the file and get a field out */
  class libconfig_string_query libconfig_string_query( LIBCONFIG_UNIT_TEST_FILE );

  SECTION( "read raw libconfig data" )
  {
    std::string query = "surfaces.flux.E";

    double e = 0;
    
    libconfig_string_query( query, e );

    REQUIRE( e == 1000 );
  }

  /* build a config module similar to how you would in a real situation and test it */
  /* Just get this to compile and get the skeleton set up. Commit that work and
     then you can stop for today */
  SECTION( "config_module" )
  {
    std::string config_path = "impurityParticleSource"

    /* constructor should create all the submodules necessary underneath itself */
    class impurity_particle_source config_module( config_path );

    /* get a value from this config module */
    int source_material_z = 0;

    config_module.get( config_module::source_material_z, source_material_z );

    REQUIRE( source_material_z == 13 )

    /* get a child module and obtain a value from it */
    auto ionization = config_module.get( config_module::ionization );

    std::string dense_grid = "";

    ionization.get( ionization::dense_grid, dense_grid );

    REQUIRE( dense_grid == "gridDensity_Ionization" );
  }

  /* for the final test, create a top-level class that includes this one and another one */

}
