#include <iostream>
#include <thrust/execution_policy.h>
#include "test_utils.hpp"
#include "config_interface.h"
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "boris.h"
#include "Surfaces.h"
#include "geometryCheck.h"

TEST_CASE( "boris" )
{
  SECTION( "e cross b" )
  {
    /* timesteps */
    int nT = 1e6;
    //int nT = 1;

    /* create particles */
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, BORIS_TEST_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    /* 2nd argument is deprecated - random number stream related */
    auto particleArray = new Particles( 1, 1, cfg_geom, gitr_flags );

    /* Captain! Important equation to convert eV energy to vector velocity components */
    gitr_precision E = 14;
    gitr_precision amu = 27;
    gitr_precision vtotal = std::sqrt(2.0 * E * 1.602e-19 / amu / 1.66e-27);
    std::cout << "vtotal: " << vtotal << std::endl;

    /* set particle properties: */
    gitr_precision dt = 1.0e-6;
    particleArray->setParticleV( 0, 0, 0, 0, vtotal, 0, 0, 13, amu, 2.0, dt );
    /*
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt)
    */

    thrust::counting_iterator<std::size_t> particle_iterator_start(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(1);

    /* Captain... Just rip everything from cross field diffusion tests and convert it
       into a boris test. Compare and contrast differences between what's in gitr.cpp
       and what's in cross_field_diffusion_tests.cpp */
    int nLines = 2;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );

    REQUIRE( nSurfaces == 2 );

    /* start */
    int nHashes = 1;
    sim::Array<int> nR_closeGeom(nHashes, 0);
    sim::Array<int> nY_closeGeom(nHashes, 0);
    sim::Array<int> nZ_closeGeom(nHashes, 0);
    sim::Array<int> nHashPoints(nHashes, 0);
    sim::Array<int> n_closeGeomElements(nHashes, 0);
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
    /* end */
    /* 
       USE_ADAPTIVE_DT = 0 run with both options, should get same answer. Start with 0 
       turn off all hashing stuff
    */

    gitr_precision perpDiffusionCoeff = 0.0;
    cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                          perpDiffusionCoeff);

    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    /* required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1 */
    /* create a unified setup script */
    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    /* uniform bfield */
    br[ 0 ] = 0;
    /* large bfield in teslas gives smaller gyromotion radius */
    by[ 0 ] = 5;
    bz[ 0 ] = 0;

    /* for the uniform efield, set efield to 1000 in z just make the cross product geometry */
    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

    int nR_PreSheathEfield = 1;
    int nY_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;

    int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

    sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);

    sim::Array<gitr_precision> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);
    PSEr[ 0 ] = 0;
    PSEz[ 0 ] = -1000;
    /* y and t */
    PSEt[ 0 ] = 0;

    int nR_closeGeom_sheath = 1;
    int nY_closeGeom_sheath = 1;
    int nZ_closeGeom_sheath = 1;
    int n_closeGeomElements_sheath = 1;
    int nGeomHash_sheath = 1;
    sim::Array<gitr_precision> closeGeomGridr_sheath(nR_closeGeom_sheath),
      closeGeomGridy_sheath(nY_closeGeom_sheath),
      closeGeomGridz_sheath(nZ_closeGeom_sheath);
    sim::Array<int> closeGeom_sheath(nGeomHash_sheath);

    /* create boris operator */
    move_boris boris( particleArray, dt, boundaries.data(), nLines, nR_Bfield, nZ_Bfield,
                      bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(), by.data(),
                      nR_PreSheathEfield, 
                      nY_PreSheathEfield,
                      nZ_PreSheathEfield,
                      &preSheathEGridr.front(), &preSheathEGridy.front(),
                      &preSheathEGridz.front(), &PSEr.front(), &PSEz.front(), &PSEt.front(),
                      nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
                      n_closeGeomElements_sheath, closeGeomGridr_sheath.data(),
                      &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
                      &closeGeom_sheath.front(), gitr_flags );

    geometry_check geometry_check0(
        particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
        nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
        n_closeGeomElements.data(), &closeGeomGridr.front(),
        &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
        nEdist, E0dist, Edist, nAdist, A0dist, Adist);

    /* get particle xyz before */
    /* time loop */
    std::cout << "Captain! num particles: " << particleArray->nParticles << std::endl;
    std::cout << "Captain! Before: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;

    for (int tt = 0; tt < nT; tt++)
    {

      thrust::for_each( thrust::device,
                        particle_iterator_start,
                        particle_iterator_end,
                        boris );

      thrust::for_each(thrust::device,
                       particle_iterator_start, 
                       particle_iterator_end, 
                       geometry_check0 );
    }

    std::cout << "Captain! After: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;
    /* get particle xyz after */

    /* what should it be analytically? */

    /* charge coulombs, e in V/m, b in Teslas */
    /* q * e / b^2 */
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
