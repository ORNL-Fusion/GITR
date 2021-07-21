#include <iostream>
#include <thrust/execution_policy.h>
#include "test_utils.hpp"
#include "config_interface.h"
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"

TEST_CASE( "boris" )
{
  SECTION( "e cross b" )
  {
    /* timesteps */
    int nT = 1e4;

    /* create particles */
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, CROSS_FIELD_GEOM_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    auto particleArray = new Particles( 1, 1, cfg_geom, gitr_flags );

    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(1);

    /* Captain... Just rip everything from cross field diffusion tests and convert it
       into a boris test. Compare and contrast differences between what's in gitr.cpp
       and what's in cross_field_diffusion_tests.cpp */
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );

    REQUIRE( nSurfaces == 2 );

    auto particleArray = new Particles( nP, 1, cfg_geom, gitr_flags );

    int nEdist = 1;
    gitr_precision E0dist = 0.0;
    gitr_precision Edist = 0.0;
    int nAdist = 1;
    gitr_precision A0dist = 0.0;
    gitr_precision Adist = 0.0;

    auto surfaces = new Surfaces(nSurfaces, nEdist, nAdist);
    surfaces->setSurface(nEdist, E0dist, Edist, nAdist, A0dist, Adist);
    sim::Array<gitr_precision> closeGeomGridr(1),
      closeGeomGridy(1), closeGeomGridz(1);
    sim::Array<int> closeGeom(1, 0);

    gitr_precision perpDiffusionCoeff = 0.0;
    cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                          perpDiffusionCoeff);

    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
  }
}
/*

necessary components:

1. particles need:

everything in impurityParticleSource.initialConditions gitrInput.cfg from protoMPex example.
but nothing other than that module.


2. boris

3. placeholder geometry


flags to set/unset:

one particle is all that should be needed for this test



boris captures drifts - which is caused by rotating particles in fields that make them rotate?
and
normal Lorentz motion - from force fields acting on it

but not

diffusion - particle movement from collisions without explicitly modeling the collision event

final answer of this test: will be drift velocity


this test essentially sets up a uniform E and B field, and tests that the boris operator
correctly simulates the drift motion of the particles. It can be a one particle test.
Simply measure the start and end position, divide by time, get drift velocity, plug into the
equations and see if it equals the analytical solution

geometry related, checking algorithms. Respecting geometry stuff.


Priorities:

coupled simulations with SOLPS-ITER,
modularize configuration and geometry components ---> how?
make it compile for mac and linux and pass tests
finish unit tests
come up with new ones
finish changing all the build time options into runtime options
What options can I take out? Which ones are defunct?

Alyssa's stuff:

1. get your version of the repo sorted out and incorporated into the new one
2. make sure GITR_legacy_python works without issue with GITR.

2015 had an IPS fusion milestone etc. Is IPS that hard to use? New solution needed.
Let's find some other alternatives for coupled simulation alternatives.

How to automate it? Talk to Wael Elwasif.

boris efield and geometry and cylindrical options.
efield test is done in matlab
















*/
