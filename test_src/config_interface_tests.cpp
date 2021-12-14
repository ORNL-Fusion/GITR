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

    /* create testing dummy */
    class testing_dummy final : public config_module_base
    {
      public:

      enum : int
      {
        setting_doubles,
        setting_string,
        nonexistent_setting,
        unregistered
      };

      testing_dummy( class libconfig_string_query const &query,
                    std::string module_path = "testing_dummy" )
        :
        config_module_base( query, module_path )
      {
        lookup[ testing_dummy::setting_doubles ] = "setting_doubles";
        lookup[ testing_dummy::setting_string ] = "setting_string";
        lookup[ testing_dummy::nonexistent_setting ] = "not_there_string_key";

        /* lookup[ testing_dummy::unregistered ] is left unset to test exception triggering */
      }
    };

    /* trigger exceptions to test for error conditions */
    /* look up an unregistered config string: */
    auto testing_dummy = std::make_shared< class testing_dummy >( query );

    /* test vectors */
    auto doubles = testing_dummy->get< std::vector< double > >( testing_dummy::setting_doubles );

    REQUIRE( doubles == std::vector< double >{ 1.5, 2.5, 3.5, 4.5, 5.5 } );

    /* test looking up unregistered key */
    try
    {
      auto trash = testing_dummy->get( testing_dummy::unregistered );
    }

    catch( class unregistered_config_mapping const &exception )
    {
      std::string const error_message{ exception.what() };

      std::cout << error_message << exception.get_key() << std::endl;
    }

    /* test looking up an unregistered key */
    /*
    REQUIRE_THROWS_MATCHES( testing_dummy->get( testing_dummy::unregistered ),
                            unregistered_config_mapping,
                            Catch::Matchers::Message( 
                            unregistered_config_mapping::get_message() ) ); 
    */

    /* test lookup up registered key with non-existent config value */
    bool caught = false;
    std::string key;

    try
    {
      int trash = testing_dummy->get< int >( testing_dummy::nonexistent_setting );
    }

    catch( class invalid_key const &exception )
    {
      caught = true;

      std::string const error_message{ exception.what() };

      std::string const error_key{ exception.get_key() };

      REQUIRE( error_key == "testing_dummy.not_there_string_key" );

      REQUIRE( exception.what() == exception.get_message() );

      std::cout << error_message << error_key << std::endl;
    }

    REQUIRE( caught ); 

    /* Captain! Make a PR for all of this so it is visible on the pull requests in GITR */
    /* test lookup on vector type vs scalar type */
    /* test lookup on wrong data type? Will that actually fail or just be a mismatch? */
    /* test all this with multiple arrows, and also from fewer levels down */
    /*

    1. Practice problems
    2. Vim
    3. Atul
    4. Complete config interface tests
    5. Fix and test operators in cross field diffusion tests, and document what you fixed
    6. Install GITR on this IDA machine

    move it to "base" and "derived" as well
    what config interface tests are still needed?

    errors we need to trigger:

    1. invalid key
    2. lookup failed - try mismatched types maybe?
    3. unregistered config mapping

    */

    /* try opening a bad file too. After this, be done with GITR this week please */

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












