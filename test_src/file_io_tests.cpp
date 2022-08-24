//#include "file_io.hpp"
#include <iostream>
#include <libconfig.h++>
#include <stdio.h>
#include "utils.h"
#include "catch2/catch_all.hpp"
#include "Particles.h"
#include "geometryCheck.h"
#include "test_data_filepath.hpp"
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
TEST_CASE("Factorials are computed", "[factorial]") {

  int const flux_ea = 0;
  int const surface_model = 0;
  int const use_3d_geom = 0;
  int const cylsymm = 0;
  int const geom_hash = 0;
  int const surface_potential = 0;

  SECTION("int")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test_data/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, DATA_TYPES_TEST_FILE );
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.thing";
    int this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    REQUIRE(this_thing == 1);
  }
  SECTION("gitr_precision")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test_data/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, DATA_TYPES_TEST_FILE );
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.float";
    gitr_precision this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    gitr_precision tol = 1e-3;
    REQUIRE_THAT(this_thing,
                     Catch::Matchers::WithinAbs(0.12345, tol));
  }
  SECTION("double")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test_data/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, DATA_TYPES_TEST_FILE );
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.float";
    double this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    double tol = 1e-3;
    REQUIRE_THAT(this_thing,
                     Catch::Matchers::WithinAbs(0.12345, tol));
  }
  SECTION("string")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test_data/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, DATA_TYPES_TEST_FILE );
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.filename";
    std::string this_thing;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    REQUIRE(this_thing == "netcdf_file_py.nc");
  }
  SECTION("geom test")
  {
  libconfig::Config cfg,cfg_geom;
  
    importLibConfig(cfg_geom, GEOM_TEST_FILE );
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  std::cout << "Start of geometry import" << std::endl;
  int nLines = 1;
  int nSurfaces = 0;
    try {
      libconfig::Setting &geom = cfg_geom.lookup("geom");
      std::cout << "Got geom setting" << std::endl;
      nLines = geom["x1"].getLength();
      std::cout << "Just read nLines " << nLines << std::endl;
      std::cout << "Number of Geometric Objects To Load: " << nLines
                << std::endl;
    } catch (const libconfig::SettingNotFoundException &nfex) {
      std::cerr << "No 'geom' setting in configuration file." << std::endl;
    }
  
    sim::Array<Boundary> boundaries(nLines + 1, Boundary());
    nSurfaces = importGeometry(cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );
    std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces
              << std::endl;
    int nParticles = 1;
  auto gitr_flags = new Flags(cfg);
      auto particleArray = new Particles(nParticles,1,cfg,gitr_flags);
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
  gitr_precision dt = 1.0e9;
  geometry_check geometry_check0(
      particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
      nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
      n_closeGeomElements.data(), &closeGeomGridr.front(),
      &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
      nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
      geom_hash,
      use_3d_geom,
      cylsymm );
  }
  
  SECTION("geom test - flat line")
  {
  libconfig::Config cfg, cfg_geom;
  
    importLibConfig(cfg_geom, FLAT_LINE_TEST_FILE );
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  std::cout << "Start of geometry import" << std::endl;
  int nLines = 1;
  int nSurfaces = 0;
    try {
      libconfig::Setting &geom = cfg_geom.lookup("geom");
      std::cout << "Got geom setting" << std::endl;
      nLines = geom["x1"].getLength();
      std::cout << "Just read nLines " << nLines << std::endl;
      std::cout << "Number of Geometric Objects To Load: " << nLines
                << std::endl;
    } catch (const libconfig::SettingNotFoundException &nfex) {
      std::cerr << "No 'geom' setting in configuration file." << std::endl;
    }
  
    sim::Array<Boundary> boundaries(nLines + 1, Boundary());
    nSurfaces = importGeometry(cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );
    std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces
              << std::endl;
    int nParticles = 1;
  auto gitr_flags = new Flags(cfg);
      auto particleArray = new Particles(nParticles,1,cfg,gitr_flags);
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
  gitr_precision dt = 1.0e9;
  geometry_check geometry_check0(
      particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
      nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
      n_closeGeomElements.data(), &closeGeomGridr.front(),
      &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
      nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
      geom_hash,
      use_3d_geom,
      cylsymm );

    particleArray->x[0] = 0.5;
    particleArray->xprevious[0] = 0.5;
    particleArray->z[0] = 0.5;
    particleArray->zprevious[0] = -0.5;

    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    std::cout << "particle hitWall " << particleArray->hitWall[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);
    
    particleArray->x[0] = 1.0;
    particleArray->xprevious[0] = 0.0;
    particleArray->z[0] = 0.5;
    particleArray->zprevious[0] = -0.15;
    particleArray->hitWall[0] = 0;
    
    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);
    
    particleArray->x[0] = 0.1;
    particleArray->xprevious[0] = 1.0;
    particleArray->y[0] = 0.0;
    particleArray->yprevious[0] = 0.0;
    particleArray->z[0] = 0.5;
    particleArray->zprevious[0] = -0.5;
    particleArray->hitWall[0] = 0;
    
    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);

  }
  
  SECTION("geom test - positive slope")
  {
  libconfig::Config cfg,cfg_geom;
  
    importLibConfig(cfg_geom, POSITIVE_SLOPE_TEST_FILE );
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  std::cout << "Start of geometry import" << std::endl;
  int nLines = 1;
  int nSurfaces = 0;
    try {
      libconfig::Setting &geom = cfg_geom.lookup("geom");
      std::cout << "Got geom setting" << std::endl;
      nLines = geom["x1"].getLength();
      std::cout << "Just read nLines " << nLines << std::endl;
      std::cout << "Number of Geometric Objects To Load: " << nLines
                << std::endl;
    } catch (const libconfig::SettingNotFoundException &nfex) {
      std::cerr << "No 'geom' setting in configuration file." << std::endl;
    }
  
    sim::Array<Boundary> boundaries(nLines + 1, Boundary());
    nSurfaces = importGeometry(cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );
    std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces
              << std::endl;
    int nParticles = 1;
  auto gitr_flags = new Flags(cfg);
      auto particleArray = new Particles(nParticles,1,cfg,gitr_flags);
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
  gitr_precision dt = 1.0e9;
  geometry_check geometry_check0(
      particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
      nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
      n_closeGeomElements.data(), &closeGeomGridr.front(),
      &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
      nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
      geom_hash,
      use_3d_geom,
      cylsymm );

    particleArray->x[0] = 0.5;
    particleArray->xprevious[0] = 0.5;
    particleArray->z[0] = 1.0;
    particleArray->zprevious[0] = -0.5;

    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);
    
    particleArray->x[0] = 0.8;
    particleArray->xprevious[0] = 0.2;
    particleArray->z[0] = 0.55;
    particleArray->zprevious[0] = 0.45;
    particleArray->hitWall[0] = 0;
    
    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);
    
    particleArray->x[0] = 0.9;
    particleArray->xprevious[0] = 0.1;
    particleArray->y[0] = 0.0;
    particleArray->yprevious[0] = 0.0;
    particleArray->z[0] = -0.5;
    particleArray->zprevious[0] = 0.5;
    particleArray->hitWall[0] = 0;
    
    geometry_check0(0);
    std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;
    REQUIRE(particleArray->hitWall[0] == 1);

  }
}
