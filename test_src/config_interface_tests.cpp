#include <iostream>
#include "test_utils.hpp"
#include "config_interface.h"

/* contains filepath of unit testing files */
#include "test_data_filepath.hpp"

TEST_CASE( "Simulation Configuration" )
{
  /* section test config modules - test instantiations from the example config file */
  SECTION( "Instantiate configuration" )
  {
    INFO( "Ahoy, Captain! File location: " << LIBCONFIG_UNIT_TEST_FILE )
    REQUIRE( 1 == 0 );
  }

  /* section libconfig data reader */

  /* section for command line parser class */

  /* This is the final test - create a config class out of an example file and read in the
     values */
}
