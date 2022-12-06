#include <iostream>
#include "catch2/catch_all.hpp"
#include "crossFieldDiffusion.h"
#include "test_data_filepath.hpp"
#include "Particles.h"
#include "geometryCheck.h"
#include "spectroscopy.h"
#include "curandInitialize.h"
#include <thrust/execution_policy.h>
#include "config_interface.h"

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T epsilon, T margin)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) 
  {
    
    bool margin_check = (a[i] != Catch::Approx(b[i]).margin(margin));
    bool epsilon_check = (a[i] != Catch::Approx(b[i]).epsilon(epsilon));
    
    if (margin_check && epsilon_check)
    {
      
      std::cout << "margin epsilon " <<
        margin_check << " " << epsilon_check << std::endl; 
      std::cout << "Element " << i << 
        " " << a[i] << " Should == " << b[i] << std::endl;
      
      return false;
    }
  }
  
  return true;
}

TEST_CASE( "cross-field diffusion operator" )
{
  /*

    Note: for this test to work, the project must be compiled with 

    GITR_USE_PERP_DIFFUSION=1
    USE_CYLSYMM=0

  */

  int const flux_ea = 1;
  int const surface_model = 1;
  int const bfield_interp = 0;
  int const use_3d_geom = 0;
  int const geom_hash = 0;
  int const spectroscopy = 2;
  int const surface_potential = 0;

  SECTION( "straight" )
  {
    int const perp_diffusion = 1;
    int const cylsymm = 0;

    /* timesteps */
    int nT = 500000;

    gitr_precision dt = 2.0e-7;

    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, CROSS_FIELD_GEOM_FILE);

    /* Captain! Instantiate a "use" module here? Get the config values following... */
    class libconfig_string_query query( CROSS_FIELD_GEOM_FILE );

    /* Captain! This will replace "flags" */
    auto use = std::make_shared< class use >( query );

    auto gitr_flags = new Flags( cfg_geom );

    int nLines = 0;

    libconfig::Setting &geom = cfg_geom.lookup( "geom" );

    nLines = geom["x1"].getLength();

    REQUIRE( nLines == 2 );

    /* particles */
    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    int nP = impurity[ "nP" ];

    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    /* Ahoy, Captain! Function call, drop in, why don't ye?! */
    int nSurfaces = importGeometry( cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );

    REQUIRE( nSurfaces == 2 );

    /* Ahoy, Captain! Function call, drop in, why don't ye?! */
    auto particleArray = new Particles( nP, 1, cfg_geom, gitr_flags );
    std::cout << "p " << particleArray->charge[0] << " x " << particleArray->xprevious[0]<<std::endl;

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

    /* Captain! Tons of build_time options in geometry_check */
    geometry_check geometry_check0(
        particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
        nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
        n_closeGeomElements.data(), &closeGeomGridr.front(),
        &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
        nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
        geom_hash,
        use_3d_geom,
        cylsymm );


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

    sim::Array<gitr_precision> gridX_midpoints(net_nX-1);
    sim::Array<gitr_precision> net_Bins((nBins + 1) * net_nX * net_nZ, 0.0);
//net_Bins needs adjusting to be in the middle of the grids
    for (int i = 0; i < net_nX; i++) {
      gridX_bins[i] = netX0 + 1.0 / (net_nX - 1) * i * (netX1 - netX0);
      std::cout << i << " " << gridX_bins[i] << std::endl;
    }

    double dx = gridX_bins[1] - gridX_bins[0];
    for (int i = 0; i < net_nX-1; i++) {
      gridX_midpoints[i] = gridX_bins[i] + 0.5*dx;
      std::cout << i << " " << gridX_midpoints[i] << std::endl;
    }

    for (int i = 0; i < net_nZ; i++) {
      gridZ_bins[i] = netZ0 + i * 1.0 / (net_nZ - 1) * (netZ1 - netZ0);
    }
    spec_bin spec_bin0(gitr_flags,particleArray, nBins, net_nX, net_nY, net_nZ,
                     &gridX_bins.front(), &gridY_bins.front(),
                     &gridZ_bins.front(), &net_Bins.front(), dt, cylsymm, spectroscopy );

    #if USE_CUDA > 0
    typedef curandState rand_type;
    #else
    typedef std::mt19937 rand_type;
    #endif

    sim::Array<rand_type> state1(nP);
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nP);
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<rand_type>(&state1.front(), true));

  gitr_precision perpDiffusionCoeff = 0.0;
  cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                        perpDiffusionCoeff);

  int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

  sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

  gitr_precision zero = 0;
  std::string empty = "";
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  importVectorField(cfg_geom, "", bfield_interp, bfieldCfg, nR_Bfield,
      0, nZ_Bfield, bfieldGridr.front(),
      zero, bfieldGridz.front(), br.front(),
      by.front(), bz.front(), empty );

  crossFieldDiffusion crossFieldDiffusion0( gitr_flags,
      particleArray, dt, &state1.front(), perpDiffusionCoeff, nR_Bfield,
      nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(),
      &bz.front(), &by.front(), perp_diffusion, cylsymm );
    
  // half-side length
    double s = 0.2;
    std::vector<gitr_precision> gold(net_nX-1,0.0);
    for( int x_bin = 0; x_bin < net_nX-1; ++x_bin )
    {
      std::cout << "xbin " << x_bin << std::endl;
      gold[x_bin] = 0.5/perpDiffusionCoeff*(s-gridX_midpoints[x_bin]);
    }

    for (int tt = 0; tt < nT; tt++)
    {
      if(tt%(nT/10) == 0) std::cout << 100.0f*tt/nT << " % done" << std::endl;
      /* call spec_bin */
      thrust::for_each(thrust::device,
                      particle_iterator0, 
                       particle_iterator_end, 
                       spec_bin0 );

      /* call diffusion */
      thrust::for_each(thrust::device,
                       particle_iterator0, 
                       particle_iterator_end, 
                       crossFieldDiffusion0 );

      /* call geometry check */
      thrust::for_each(thrust::device,
                       particle_iterator0, 
                       particle_iterator_end, 
                       geometry_check0 );
    }

    /* examine histogram */
    /* output x */
    /* analytical equation we expect */
    std::vector<double> density(net_nX-1,0.0);
    double sum = 0;
    double rms_error = 0;
    for( int x_bin = 0; x_bin < net_nX-1; ++x_bin )
    {
      /* sum over z? */
      for( int z_bin = 0; z_bin < net_nZ; ++z_bin )
      {
        sum += net_Bins[ nBins * net_nX * net_nZ +
                         z_bin * net_nX + x_bin ];
      }

      //density[x_bin] = sum;
      density[x_bin] = sum*dt/nP/dx;
      std::cout << "dens: " << density[x_bin] << " gold " << gold[x_bin] << std::endl;
      rms_error = rms_error + (density[x_bin] - gold[x_bin])*(density[x_bin] - gold[x_bin]);
      sum = 0.0;
    }

    rms_error = std::sqrt(rms_error/(net_nX-1));
    std::cout << "rms_error " << rms_error << std::endl;

    int not_hit = 0;    
    for (int i=0; i< nP; i++ )
    {
      if(particleArray->hitWall[i]<1.0)
      {
        not_hit = not_hit+1;
      }
    //    std::cout << "p " << particleArray->charge[i] << " x " << particleArray->xprevious[i]<<std::endl;
    }   

    std::cout << "number of particles not finished " << not_hit << std::endl;
    gitr_precision margin = 0.1;
    gitr_precision epsilon = 0.05;
    REQUIRE(rms_error<6.0e-4);
  }

  /*

    Note: for this test to work, the project must be compiled with 

    GITR_USE_PERP_DIFFUSION=2
    USE_CYLSYMM=1

    USE3DTETGEOM=0
    USE_SURFACE_POTENTIAL=0
    PARTICLE_SOURCE_FILE=0
    USE_ADAPTIVE_DT=0

  */
  SECTION( "curved" )
  {
    int const perp_diffusion = 2;
    /* timesteps */
    int nT = 200000;
    int const cylsymm = 1;

    gitr_precision dt = 1.0e-7;

    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, CROSS_FIELD_GEOM_FILE_1 );

    auto gitr_flags = new Flags( cfg_geom );

    int nLines = 0;

    libconfig::Setting &geom = cfg_geom.lookup( "geom" );

    nLines = geom["x1"].getLength();

    REQUIRE( nLines == 2 );

    /* particles */
    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    int nP = impurity[ "nP" ];

    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );

    REQUIRE( nSurfaces == 2 );

    auto particleArray = new Particles( nP, 1, cfg_geom, gitr_flags );
    for(int i=0;i<nP;i++)
    {
      particleArray->setParticleV(i,0.2,0.0,0.0, 10, 10, 10, 4, 5, 1.0,dt);
    }
    std::cout << "p " << particleArray->charge[0] << " x " << particleArray->xprevious[0]<<std::endl;

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
        nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
        geom_hash,
        use_3d_geom,
        cylsymm );


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

    sim::Array<gitr_precision> gridX_midpoints(net_nX-1);
    sim::Array<gitr_precision> net_Bins((nBins + 1) * net_nX * net_nZ, 0.0);
//net_Bins needs adjusting to be in the middle of the grids
    for (int i = 0; i < net_nX; i++) {
      gridX_bins[i] = netX0 + 1.0 / (net_nX - 1) * i * (netX1 - netX0);
      std::cout << i << " " << gridX_bins[i] << std::endl;
    }

    double dx = gridX_bins[1] - gridX_bins[0];
    for (int i = 0; i < net_nX-1; i++) {
      gridX_midpoints[i] = gridX_bins[i] + 0.5*dx;
      std::cout << i << " " << gridX_midpoints[i] << std::endl;
    }

    for (int i = 0; i < net_nZ; i++) {
      gridZ_bins[i] = netZ0 + i * 1.0 / (net_nZ - 1) * (netZ1 - netZ0);
    }
    spec_bin spec_bin0(gitr_flags,particleArray, nBins, net_nX, net_nY, net_nZ,
                     &gridX_bins.front(), &gridY_bins.front(),
                     &gridZ_bins.front(), &net_Bins.front(), dt, cylsymm, spectroscopy );

    #if USE_CUDA > 0
    typedef curandState rand_type;
    #else
    typedef std::mt19937 rand_type;
    #endif

    sim::Array<rand_type> state1(nP);
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nP);
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<rand_type>(&state1.front(), true));

  gitr_precision perpDiffusionCoeff = 0.0;
  cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                        perpDiffusionCoeff);

  int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

  sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

  gitr_precision zero = 0;
  std::string empty = "";
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  importVectorField(cfg_geom, "", bfield_interp, bfieldCfg, nR_Bfield,
      0, nZ_Bfield, bfieldGridr.front(),
      zero, bfieldGridz.front(), br.front(),
      by.front(), bz.front(), empty );

  crossFieldDiffusion crossFieldDiffusion0( gitr_flags,
      particleArray, dt, &state1.front(), perpDiffusionCoeff, nR_Bfield,
      nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(),
      &bz.front(), &by.front(), perp_diffusion, cylsymm );
    
  // half-side length
    double s = 0.2;
    std::vector<gitr_precision> gold(net_nX-1,0.0);
    for( int x_bin = 0; x_bin < net_nX-1; ++x_bin )
    {
      std::cout << "xbin " << x_bin << std::endl;
      //flux*a/D*log(b/r)
      gold[x_bin] = 0.5*0.2/perpDiffusionCoeff*std::log(0.3/gridX_midpoints[x_bin]);
    }

    for (int tt = 0; tt < nT; tt++)
    {
      if(tt%(nT/10) == 0) std::cout << 100.0f*tt/nT << " % done" << std::endl;
      /* call spec_bin */
      thrust::for_each(thrust::device,
                      particle_iterator0, 
                       particle_iterator_end, 
                       spec_bin0 );

      /* call diffusion */
      thrust::for_each(thrust::device,
                       particle_iterator0, 
                       particle_iterator_end, 
                       crossFieldDiffusion0 );

      /* call geometry check */
      thrust::for_each(thrust::device,
                       particle_iterator0, 
                       particle_iterator_end, 
                       geometry_check0 );
    }

    /* examine histogram */
    /* output x */
    /* analytical equation we expect */
    std::vector<double> density(net_nX-1,0.0);
    double sum = 0;
    double rms_error = 0;
    for( int x_bin = 0; x_bin < net_nX-1; ++x_bin )
    {
      /* sum over z? */
      for( int z_bin = 0; z_bin < net_nZ; ++z_bin )
      {
        sum += net_Bins[ nBins * net_nX * net_nZ +
                         z_bin * net_nX + x_bin ];
      }
      //density[x_bin] = sum;
      density[x_bin] = 0.5/3.1415*sum*dt/nP/gridX_midpoints[x_bin]/dx;
      std::cout << "dens: " << density[x_bin] << " gold " << gold[x_bin] << std::endl;
      rms_error = rms_error + (density[x_bin] - gold[x_bin])*(density[x_bin] - gold[x_bin]);
      sum = 0.0;
    }
    rms_error = std::sqrt(rms_error/(net_nX-1));
    std::cout << "rms_error " << rms_error << std::endl;

int not_hit = 0;    
for (int i=0; i< nP; i++ )
{
  if(particleArray->hitWall[i]<1.0)
  {
    not_hit = not_hit+1;
  }
//    std::cout << "p " << particleArray->charge[i] << " x " << particleArray->xprevious[i]<<std::endl;
}   
    std::cout << "number of particles not finished " << not_hit << std::endl;
    gitr_precision margin = 0.1;
    gitr_precision epsilon = 0.05;
    REQUIRE(rms_error<5.0e-4);
  }
}
