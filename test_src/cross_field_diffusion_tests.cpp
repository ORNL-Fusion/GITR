#include <iostream>
#include "test_utils.hpp"
#include "crossFieldDiffusion.h"
#include "test_data_filepath.hpp"
#include "Particles.h"
#include "geometryCheck.h"

TEST_CASE( "cross-field diffusion operator" )
{
  SECTION( "Ahoy, Captain!" )
  {
    /* timesteps */
    int nT = 1e4;

    /* particles */
    int nP = 1;

    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, CROSS_FIELD_GEOM_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    int nLines = 0;

    libconfig::Setting &geom = cfg_geom.lookup( "geom" );

    nLines = geom["x1"].getLength();

    REQUIRE( nLines == 2 );

    /* Correct flags for this unit test */
    /*
    USE3DTETGEOM=0
    USE_SURFACE_POTENTIAL=0
    PARTICLE_SOURCE_FILE=0
    */

    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );

    REQUIRE( nSurfaces == 2 );

    auto particleArray = new Particles( nP, 1, cfg_geom, gitr_flags );
  }
}
