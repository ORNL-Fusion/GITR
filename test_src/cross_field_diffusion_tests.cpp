#include <iostream>
#include "test_utils.hpp"
#include "crossFieldDiffusion.h"
#include "test_data_filepath.hpp"
#include "Particles.h"
#include "geometryCheck.h"
#include "spectroscopy.h"
#include "curandInitialize.h"
#include <thrust/execution_policy.h>

/* is there an option for analytic bfield? No. */
TEST_CASE( "cross-field diffusion operator" )
{
  SECTION( "Ahoy, Captain!" )
  {
    /* timesteps */
    int nT = 1e4;

    gitr_precision dt = 1.0e-6;

    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, CROSS_FIELD_GEOM_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    int nLines = 0;

    libconfig::Setting &geom = cfg_geom.lookup( "geom" );

    nLines = geom["x1"].getLength();

    REQUIRE( nLines == 2 );

    /* particles */
    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    int nP = impurity[ "nP" ];

    /* Correct flags for this unit test */
    /*
    USE3DTETGEOM=0
    USE_SURFACE_POTENTIAL=0
    PARTICLE_SOURCE_FILE=0
    SPECTROSCOPY=2 ( for 2D, 1 is not an option )
    USEPERPDIFFUSION=1 (PARDIFFUSION is deprecated )
    USE_ADAPTIVE_DT=0
    */

    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );

    REQUIRE( nSurfaces == 2 );

    auto particleArray = new Particles( nP, 1, cfg_geom, gitr_flags );

    /* dummy variables for hashing */
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

    geometry_check geometry_check0(
        particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
        nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
        n_closeGeomElements.data(), &closeGeomGridr.front(),
        &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
        nEdist, E0dist, Edist, nAdist, A0dist, Adist);


    /* data collection variables */
    gitr_precision netX0 = 0.0;
    gitr_precision netX1 = 0.0, netY0 = 0.0, netY1 = 0.0, netZ0 = 0.0, netZ1 = 0.0;
    int net_nX = 0, net_nY = 0, net_nZ = 0;
    int nBins = 0;

    /* populate from diagnostics */
    if (cfg_geom.lookupValue("diagnostics.netx0", netX0) &&
        cfg_geom.lookupValue("diagnostics.netx1", netX1) &&
        cfg_geom.lookupValue("diagnostics.netz0", netZ0) &&
        cfg_geom.lookupValue("diagnostics.netz1", netZ1) &&
        cfg_geom.lookupValue("diagnostics.nX", net_nX) &&
        cfg_geom.lookupValue("diagnostics.nZ", net_nZ) &&
        cfg_geom.lookupValue("diagnostics.densityChargeBins", nBins))
    {
      std::cout << "Spectroscopy net imported" << std::endl;
    }

    sim::Array<gitr_precision> gridX_bins(net_nX), gridY_bins(1), gridZ_bins(net_nZ);

    sim::Array<double> net_Bins((nBins + 1) * net_nX * net_nZ, 0.0);

    for (int i = 0; i < net_nX; i++) {
      gridX_bins[i] = netX0 + 1.0 / (net_nX - 1) * i * (netX1 - netX0);
    }

    for (int i = 0; i < net_nZ; i++) {
      gridZ_bins[i] = netZ0 + i * 1.0 / (net_nZ - 1) * (netZ1 - netZ0);
    }
    spec_bin spec_bin0(gitr_flags,particleArray, nBins, net_nX, net_nY, net_nZ,
                     &gridX_bins.front(), &gridY_bins.front(),
                     &gridZ_bins.front(), &net_Bins.front(), dt);

    typedef std::mt19937 rand_type;

    sim::Array<rand_type> state1(nP);
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nP);
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<>(&state1.front(), 0));

  gitr_precision perpDiffusionCoeff = 0.0;
  cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                        perpDiffusionCoeff);

  int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

  sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

  double zero = 0;
  std::string empty = "";
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  importVectorField(cfg_geom, "", BFIELD_INTERP, bfieldCfg, nR_Bfield,
      0, nZ_Bfield, bfieldGridr.front(),
      zero, bfieldGridz.front(), br.front(),
      by.front(), bz.front(), empty );

  crossFieldDiffusion crossFieldDiffusion0(gitr_flags,
      particleArray, dt, &state1.front(), perpDiffusionCoeff, nR_Bfield,
      nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(),
      &bz.front(), &by.front());

    for (int tt = 0; tt < nT; tt++)
    {
      /* call spec_bin */
      thrust::for_each(thrust::host,
                      particle_iterator0, 
                       particle_iterator_end, 
                       spec_bin0 );

      /* call diffusion */
      thrust::for_each(thrust::host,
                       particle_iterator0, 
                       particle_iterator_end, 
                       crossFieldDiffusion0 );

      /* call geometry check */
      thrust::for_each(thrust::host,
                       particle_iterator0, 
                       particle_iterator_end, 
                       geometry_check0 );
    }

    /* examine histogram */
    /* output x */
    /* analytical equation we expect */
    for( int x_bin = 0; x_bin < net_nX; ++x_bin )
    {
      /* sum over z? */
      double sum = 0;
      for( int z_bin = 0; z_bin < net_nZ; ++z_bin )
      {
        sum += net_Bins[ nBins * net_nX * net_nZ +
                         z_bin * net_nX + x_bin ];
      }
      std::cout << "sum: " << sum << std::endl;
    }
  }
}
