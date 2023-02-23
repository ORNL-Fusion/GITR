#include "file_io_data_broker.h"
#include <iostream>
#include <libconfig.h++>
#include <stdio.h>
#include "utils.h"
#include "Particles.h"
#include "geometryCheck.h"
#include "test_data_filepath.hpp"
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

file_io_data_broker::file_io_data_broker()
{
}

bool file_io_data_broker::run_0()
{
  int const flux_ea = 0;
  int const surface_model = 0;
  int const use_3d_geom = 0;
  int const cylsymm = 0;
  int const geom_hash = 0;
  int const surface_potential = 0;
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

  bool pass = ( particleArray->hitWall[0] == 1 );

  particleArray->x[0] = 1.0;
  particleArray->xprevious[0] = 0.0;
  particleArray->z[0] = 0.5;
  particleArray->zprevious[0] = -0.15;
  particleArray->hitWall[0] = 0;

  geometry_check0(0);
  std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;

  pass = pass && ( particleArray->hitWall[0] == 1 );

  particleArray->x[0] = 0.1;
  particleArray->xprevious[0] = 1.0;
  particleArray->y[0] = 0.0;
  particleArray->yprevious[0] = 0.0;
  particleArray->z[0] = 0.5;
  particleArray->zprevious[0] = -0.5;
  particleArray->hitWall[0] = 0;

  geometry_check0(0);
  std::cout << "particle x and z " << particleArray->x[0] << " " <<  particleArray->z[0] << std::endl;

  pass = pass && ( particleArray->hitWall[0] == 1 );

  return pass;
}
