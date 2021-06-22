#include <iostream>
#include "test_utils.hpp"
#include "crossFieldDiffusion.h"
#include "test_data_filepath.hpp"
#include "Particles.h"
#include "geometryCheck.h"

/* Captain! Split this into a cpp and header */
TEST_CASE( "cross-field diffusion operator" )
{
  /*
  what do all these have in common? 

  particles:
  a state
  a begin and end iterator - initialize the state
  */



  SECTION( "Ahoy, Captain!" )
  {
    /* timesteps */
    int nT = 1e4;
    /* particles */
    int nP = 1e4;

    /* create the geometry */

    /* instantiate particles */
    /* instantiate geometryCheck */
    /* instantiate */
    REQUIRE( false );
  }

  /**/
}
