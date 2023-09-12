//------------------------------------------------------------------------------
// GITR: gitr.cpp
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - August 2023 by Diaw
//     - Mass and charge are arrays now
//     - Remove unused code blocks
//     - Add new sheath model and according flags
//     - Particle source generation and reading are now in separate files: getParticleData.h
//     - Ionization and recombination reading and processing are now in separate files: processIonizationRecombination.h
//
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------


#include "Boundary.h"
//#include "Fields.h"
#include "Particles.h"
#include "Surfaces.h"
#include "array.h"
#include "boris.h"
#include "boundaryInit.h"
#include "coulombCollisions.h"
#include "crossFieldDiffusion.h"
#include "curandInitialize.h"
#include "fieldLineTrace.h"
#include "geometryCheck.h"
#include "hashGeom.h"
#include "hashGeomSheath.h"
#include "history.h"
#include "interp2d.hpp"
#include "interpolate.h"
#include "ionize.h"
#include <cmath>
#include "recombine.h"
#include "spectroscopy.h"
#include "surfaceModel.h"
#include "thermalForce.h"
#include "utils.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <netcdf.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "flags.hpp"

// AD
#include <numeric> 
#include "getParticleData.h"
//#include "surfaceReact.h"
#include "materials.h"
#ifdef __CUDACC__
#include <curand.h>
#include <curand_kernel.h>
#else
#endif

#if USE_MPI
#include <mpi.h>
#endif

// #include <omp.h>
#include </opt/homebrew/opt/libomp/include/omp.h>

#include "sortParticles.h"
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include "config_interface.h"
#include "CLI/CLI.hpp"

using namespace netCDF;


#include <vector>
#include <netcdf>




#if USE_DOUBLE
typedef double gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncDouble;
//MPI_Datatype mpi_precision = MPI_DOUBLE
#else
typedef float gitr_precision;
netCDF::NcType netcdf_precision = netCDF::ncFloat;
#endif

int main(int argc, char **argv, char **envp) {

  CLI::App app{ "!" };
  std::string file_name = "input/gitrInput.cfg";
  app.add_option( "-c", file_name, "config filepath" );

  CLI11_PARSE( app, argc, argv );

  std::cout << "file_name read from stdin: " << file_name << std::endl;

  typedef std::chrono::high_resolution_clock gitr_time;
  auto gitr_start_clock = gitr_time::now();
  class libconfig_string_query query( file_name );
  class use use( query );

  int surface_model = use.get< int >( use::surface_model );
  int flux_ea = use.get<int>(use::flux_ea);
  int spectroscopy = use.get< int >( use::spectroscopy );
  // hardcoded to 0 for now, taken out of config_interface
  int biased_surface = BIASED_SURFACE;
  int use_3d_geom = use.get< int >( use::use_3d_geom );
  int cylsymm = use.get< int >( use::cylsymm );
  int bfield_interp = use.get< int >( use::bfield_interp );
  int gradt_interp = use.get< int >( use::gradt_interp );
  int force_eval = use.get< int >( use::force_eval );
  int sort_particles = use.get< int >( use::sort );
  int use_adaptive_dt = use.get< int >( use::adaptive_dt );
  int geom_hash = use.get<int>( use::geom_hash);
  int particle_source_file = use.get< int >( use::particle_source_file );
  int particle_source_space = use.get< int >( use::particle_source_space );
  int particle_source_energy = use.get< int >( use::particle_source_energy );
  int particle_source_angle = use.get< int >( use::particle_source_angle );
  int particle_tracks = use.get< int >( use::particle_tracks );
  int presheath_interp = use.get< int >( use::presheath_interp );
  int efield_interp = use.get< int >( use::efield_interp );
  int surface_potential = use.get< int >( use::surface_potential );
  int flowv_interp = use.get<int>( use::flowv_interp );
  int density_interp = use.get<int>( use::density_interp );
  int temp_interp = use.get<int>(use::temp_interp);
  int geom_hash_sheath = use.get<int>( use::geom_hash_sheath );
  int thermal_force = use.get< int >( use::thermal_force );
  int sheath_efield = use.get< int >( use::sheath_efield );
  // hardcoded to 1 for now
  //int presheath_efield = use.get< int >( use::presheath_efield );
  int presheath_efield = 1;
  int ionization = use.get< int >( use::ionization );
  int coulomb_collisions = use.get< int >( use::coulomb_collisions );
  int perp_diffusion = use.get< int >( use::perp_diffusion );
  // hardcoded to 1 for now
  //int field_aligned_values = use.get< int >( use::field_aligned_values );
  int field_aligned_values = FIELD_ALIGNED_VALUES;
  bool fixed_seeds = bool( use.get< int >( use::fixed_seeds ) );

  // get sheath model type
  int sheath_model_type = use.get< int >( use::sheath_model_type );
  int nspecies = use.get< int >( use::nspecies );


  // Set default processes per node to 1
  int ppn = 1;

  // Set default input file string
  std::string inputFile = file_name;

#if USE_MPI > 0
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
#endif
 
  read_comand_line_args(argc,argv,ppn,inputFile);

#if USE_MPI > 0
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("\nHello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, world_rank, world_size);
#if USE_CUDA > 0
  cudaSetDevice(world_rank % ppn);
#endif
#else
  int world_rank = 0;
  int world_size = 1;
#endif
 //cudaSetDevice(1); 
  // Prepare config files for import
  libconfig::Config cfg, cfg_geom;
  cfg.setAutoConvert(true);
  cfg_geom.setAutoConvert(true);

  std::string input_path = "input/";
  
  if (world_rank == 0) {
    // Parse and read input file
    std::cout << "Open configuration file " << input_path << inputFile
              << std::endl;
    importLibConfig(cfg, inputFile);
    // Parse and read geometry file
    std::string geomFile;
    getVariable(cfg, "geometry.fileString", geomFile);
    std::cout << "Open geometry file " << input_path + geomFile << std::endl;
    importLibConfig(cfg_geom, input_path + geomFile);

    std::cout << "Successfully staged input and geometry file " << std::endl;

// check binary compatibility with input file
  }
  auto gitr_flags = new Flags(cfg);
    std::cout << "gitr flags " << gitr_flags->USE_IONIZATION << std::endl;

  // Background species info
  gitr_precision background_Z = 0.0, background_amu = 0.0;
  if (world_rank == 0) {
    getVariable(cfg, "backgroundPlasmaProfiles.Z", background_Z);
    getVariable(cfg, "backgroundPlasmaProfiles.amu", background_amu);
  }
#if USE_MPI > 0
  MPI_Bcast(&background_Z, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&background_amu, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  auto finish_clock0nc = gitr_time::now();
  typedef std::chrono::duration<gitr_precision> fsec0nc;
  fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
  //printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
  
  
  int nR_Bfield = 1, nY_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  std::string bfieldFile;
  if (world_rank == 0) {
    //std::cout << "Ahoy, Captain!" << std::endl;
    importVectorFieldNs(cfg, input_path, bfield_interp, bfieldCfg, nR_Bfield,
                        nY_Bfield, nZ_Bfield, bfieldFile);
    //std::cout << "Ahoy, Captain!" << std::endl;
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridy(nY_Bfield),
      bfieldGridz(nZ_Bfield);
  n_Bfield = nR_Bfield * nY_Bfield * nZ_Bfield;
  sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  if (world_rank == 0) {
    importVectorField(cfg, input_path, bfield_interp, bfieldCfg, nR_Bfield,
                      nY_Bfield, nZ_Bfield, bfieldGridr.front(),
                      bfieldGridy.front(), bfieldGridz.front(), br.front(),
                      by.front(), bz.front(), bfieldFile);
  }
#if USE_MPI > 0
  MPI_Bcast(br.data(), n_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(by.data(), n_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(bz.data(), n_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(bfieldGridr.data(), nR_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(bfieldGridy.data(), nY_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(bfieldGridz.data(), nZ_Bfield, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  gitr_precision Btest[3] = {0.0};
  interp2dVector(&Btest[0], 5.5, 0.0, -4.0, nR_Bfield, nZ_Bfield,
                 bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(),
                 by.data(), cylsymm );
//   std::cout << "node " << world_rank << "Bfield at 5.5 -4 " << Btest[0] << " "
//            << Btest[1] << " " << Btest[2] << std::endl;
//   std::string profiles_folder = "output/profiles";
// exit(0);

  std::string profiles_folder = "output/profiles";


  // Geometry Definition
  std::cout << "Start of geometry import" << std::endl;
  int nLines = 1;
  int nSurfaces = 0;
  if (world_rank == 0) {
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
  }


#if USE_MPI > 0
  MPI_Bcast(&nLines, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nSurfaces, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  sim::Array<Boundary> boundaries(nLines + 1, Boundary());
  if (world_rank == 0) {
    nSurfaces = importGeometry(cfg_geom, boundaries, use_3d_geom, cylsymm, surface_potential );
    std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces
              << std::endl;
  }
#if USE_MPI > 0
  MPI_Bcast(&nSurfaces, 1, MPI_INT, 0, MPI_COMM_WORLD);
  int nBoundaryMembers;
  if(use_3d_geom)
  {
    nBoundaryMembers = 41;
  }
  else
  {
    int nBoundaryMembers = 39;
  }
  int nIntMembers = 6;
  int lengths[nBoundaryMembers] = {0};
  MPI_Aint offsets[nBoundaryMembers] = {};
  MPI_Datatype types[nBoundaryMembers] = {};
  for (int i = 0; i < nBoundaryMembers; i++) {
    lengths[i] = 1;
    offsets[i] = i * 4;
    if (i < nIntMembers) {
      types[i] = MPI_INT;
    } else {
      types[i] = MPI_FLOAT;
    }
  }
  MPI_Datatype boundary_type;
  MPI_Type_create_struct(nBoundaryMembers, lengths, offsets, types,
                         &boundary_type);
  MPI_Type_commit(&boundary_type);
  MPI_Bcast(&boundaries[0], nLines + 1, boundary_type, 0, MPI_COMM_WORLD);
#endif

  gitr_precision biasPotential = 0.0;
  if(biased_surface)
  {
  if (world_rank == 0) {
    getVariable(cfg, "backgroundPlasmaProfiles.biasPotential", biasPotential);
  }
#if USE_MPI > 0
  MPI_Bcast(&biasPotential, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  // create Surface data structures
  int nEdist = 1;
  gitr_precision E0dist = 0.0;
  gitr_precision Edist = 0.0;
  int nAdist = 1;
  gitr_precision A0dist = 0.0;
  gitr_precision Adist = 0.0;
  if(flux_ea)
  {
  if (world_rank == 0) {
    getVariable(cfg, "surfaces.flux.nE", nEdist);
    getVariable(cfg, "surfaces.flux.E0", E0dist);
    getVariable(cfg, "surfaces.flux.E", Edist);

    getVariable(cfg, "surfaces.flux.nA", nAdist);
    getVariable(cfg, "surfaces.flux.A0", A0dist);
    getVariable(cfg, "surfaces.flux.A", Adist);
  }
#if USE_MPI > 0
  const int nSurfaceMembers = 18;

  int lengthsSurface[nSurfaceMembers] = {1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         nEdist * nAdist,
                                         nEdist,
                                         nAdist,
                                         nEdist * nAdist,
                                         nEdist * nAdist,
                                         nEdist * nAdist,
                                         nEdist * nAdist,
                                         nEdist * nAdist,
                                         nSurfaces * nEdist * nAdist};

  MPI_Aint offsetsSurface[nSurfaceMembers] = {
      offsetof(Surfaces, nSurfaces),
      offsetof(Surfaces, nE),
      offsetof(Surfaces, nA),
      offsetof(Surfaces, E0),
      offsetof(Surfaces, E),
      offsetof(Surfaces, A0),
      offsetof(Surfaces, A),
      offsetof(Surfaces, dE),
      offsetof(Surfaces, dA),
      offsetof(Surfaces, sumParticlesStrike),
      offsetof(Surfaces, gridE),
      offsetof(Surfaces, gridA),
      offsetof(Surfaces, sumWeightStrike),
      offsetof(Surfaces, grossDeposition),
      offsetof(Surfaces, grossErosion),
      offsetof(Surfaces, aveSputtYld),
      offsetof(Surfaces, sputtYldCount),
      offsetof(Surfaces, energyDistribution)};
  MPI_Datatype typesSurface[nSurfaceMembers] = {
      MPI_INT,   MPI_INT,   MPI_INT,   MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
      MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT,   MPI_FLOAT, MPI_FLOAT,
      MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
  MPI_Datatype surface_type;
  MPI_Type_create_struct(nSurfaceMembers, lengthsSurface, offsetsSurface,
                         typesSurface, &surface_type);
  MPI_Type_commit(&surface_type);
  MPI_Bcast(&nEdist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nAdist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&E0dist, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Edist, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&A0dist, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Adist, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  auto surfaces = new Surfaces(nSurfaces, nspecies, nEdist, nAdist);
  surfaces->setSurface(nEdist, E0dist, Edist, nAdist, A0dist, Adist);

  //#if USE_MPI > 0
  // Arrays used for reduction at end of sim
  sim::Array<gitr_precision> grossDeposition(nSurfaces * nspecies, 0.0);
  sim::Array<gitr_precision> grossErosion(nSurfaces * nspecies, 0.0);
  sim::Array<gitr_precision> sumWeightStrike(nSurfaces * nspecies, 0.0);
  sim::Array<gitr_precision> energyDistribution(nSurfaces * nspecies * nEdist * nAdist, 0.0);
  sim::Array<gitr_precision> reflDistribution(nSurfaces * nspecies * nEdist * nAdist, 0.0);
  sim::Array<gitr_precision> sputtDistribution(nSurfaces * nspecies * nEdist * nAdist, 0.0);
  sim::Array<gitr_precision> aveSputtYld(nSurfaces * nspecies, 0.0);
  sim::Array<int> sputtYldCount(nSurfaces * nspecies, 0);
  sim::Array<int> sumParticlesStrike(nSurfaces * nspecies, 0);


  int nHashes = 1;
  int nR_closeGeomTotal = 1;
  int nY_closeGeomTotal = 1;
  int nZ_closeGeomTotal = 1;
  int nHashPointsTotal = 1;
  int nGeomHash = 1;
  std::string geomHashCfg = "geometry_hash.";

  if( geom_hash == 1 )
  {
  if (world_rank == 0) {
    getVariable(cfg, geomHashCfg + "nHashes", nHashes);
  }
#if USE_MPI > 0
  MPI_Bcast(&nHashes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  sim::Array<int> nR_closeGeom(nHashes, 0);
  sim::Array<int> nY_closeGeom(nHashes, 0);
  sim::Array<int> nZ_closeGeom(nHashes, 0);
  sim::Array<int> nHashPoints(nHashes, 0);
  sim::Array<int> n_closeGeomElements(nHashes, 0);

  if( geom_hash == 1 )
  {
  if (world_rank == 0) {
    importHashNs(cfg, input_path, nHashes, "geometry_hash", nR_closeGeom.data(),
                 nY_closeGeom.data(), nZ_closeGeom.data(),
                 n_closeGeomElements.data(), nR_closeGeomTotal,
                 nY_closeGeomTotal, nZ_closeGeomTotal, nHashPoints.data(),
                 nHashPointsTotal, nGeomHash, use_3d_geom );
    std::cout << "made it here" << std::endl;
    libconfig::Setting& geomHash = cfg.lookup("geometry_hash");
     if(nHashes > 1)
    {
      for(int i=0; i<nHashes;i++)
      {
        nR_closeGeom[i] = geomHash["nR_closeGeom"][i];
        nZ_closeGeom[i] = geomHash["nZ_closeGeom"][i];
        n_closeGeomElements[i] = geomHash["n_closeGeomElements"][i];
        std::cout << "hash nr ny nz total " << n_closeGeomElements[i] << " "
        << nR_closeGeom[i]  << " " << nZ_closeGeom[i]<< std::endl;
      }
    }
     else
    {
      getVariable(cfg,geomHashCfg+"nR_closeGeom",nR_closeGeom[0]);
      getVariable(cfg,geomHashCfg+"nZ_closeGeom",nZ_closeGeom[0]);
      getVariable(cfg,geomHashCfg+"n_closeGeomElements",n_closeGeomElements[0]);
        std::cout << "else hash nr ny nz total " << n_closeGeomElements[0] << " "
        << nR_closeGeom[0]  << " " << nZ_closeGeom[0]<< std::endl;
    }
     for(int j=0;j<nHashes;j++)
    {
      nGeomHash = nGeomHash +
      nR_closeGeom[j]*nZ_closeGeom[j]*n_closeGeomElements[j];
      nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[j];
      nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[j];
    }
    if( use_3d_geom > 0 )
    {
     if(nHashes > 1)
    {
      for(int i=0; i<nHashes;i++)
      {
        nY_closeGeom[i] = geomHash["nY_closeGeom"][i];
      }
    }
     else
    {
      getVariable(cfg,geomHashCfg+"nY_closeGeom",nY_closeGeom[0]);
    }
    }
     nGeomHash = 0;
     nR_closeGeomTotal = 0;
     nY_closeGeomTotal = 0;
     nZ_closeGeomTotal = 0;
     nHashPointsTotal = 0;
     for(int j=0;j<nHashes;j++)
    {
        nHashPoints[j] =nR_closeGeom[j]*nY_closeGeom[j]*nZ_closeGeom[j];
      nHashPointsTotal = nHashPointsTotal + nHashPoints[j];
      nGeomHash = nGeomHash + nHashPoints[j]*n_closeGeomElements[j];
      nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[j];
      nY_closeGeomTotal = nY_closeGeomTotal + nY_closeGeom[j];
      nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[j];
    }
     std::cout << "hhhash nr ny nz total " << nGeomHash << " " <<
     nR_closeGeomTotal << " " << nY_closeGeomTotal << " " <<
     nZ_closeGeomTotal<< std::endl;
  }
#if USE_MPI > 0
  std::cout << " mpi broadcast hash " << std::endl;
  MPI_Bcast(&nR_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nR_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nHashPointsTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << " mpi broadcast hash finished" << std::endl;
#endif
  }


std::vector<std::string> hashFile;
if( geom_hash > 1 )
{
  if (world_rank == 0) {
    getVariable(cfg, geomHashCfg + "nHashes", nHashes);
  }
#if USE_MPI > 0
  MPI_Bcast(&nHashes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (world_rank == 0) {
	  libconfig::Setting &geomHash = cfg.lookup("geometry_hash");
    for (int i = 0; i < nHashes; i++) {
      if (nHashes > 1) {
        hashFile.push_back(geomHash["fileString"][i]);
      } else {
        hashFile.push_back("dummy");
        getVariable(cfg, geomHashCfg + "fileString", hashFile[i]);
      }
      nR_closeGeom[i] = getDimFromFile(cfg, input_path + hashFile[i],
                                       geomHashCfg, "gridNrString");
      nZ_closeGeom[i] = getDimFromFile(cfg, input_path + hashFile[i],
                                       geomHashCfg, "gridNzString");
      n_closeGeomElements[i] = getDimFromFile(
          cfg, input_path + hashFile[i], geomHashCfg, "nearestNelementsString");
      nGeomHash = nGeomHash +
                  nR_closeGeom[i] * nZ_closeGeom[i] * n_closeGeomElements[i];

      if( use_3d_geom > 0 )
      {
        nY_closeGeom[i] = getDimFromFile(cfg, input_path + hashFile[i],
                                         geomHashCfg, "gridNyString");
        nGeomHash = nGeomHash -
                    nR_closeGeom[i] * nZ_closeGeom[i] * n_closeGeomElements[i] +
                    nY_closeGeom[i] * nR_closeGeom[i] * nZ_closeGeom[i] *
                        n_closeGeomElements[i];
        nY_closeGeomTotal = nY_closeGeomTotal + nY_closeGeom[i];
      }

      nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[i];
      nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[i];
    }
  }
#if USE_MPI > 0
  std::cout << "starting mpibacast" << std::endl;
  MPI_Bcast(&nR_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements[0], nHashes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nR_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeomTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nHashPointsTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "finished mpibacast" << std::endl;
#endif
}

  std::cout << "allocating closGeomGrids " << nR_closeGeomTotal << " "
            << nY_closeGeomTotal << " " << nZ_closeGeomTotal << " " << nGeomHash
            << std::endl;
  sim::Array<gitr_precision> closeGeomGridr(nR_closeGeomTotal),
      closeGeomGridy(nY_closeGeomTotal), closeGeomGridz(nZ_closeGeomTotal);
  sim::Array<int> closeGeom(nGeomHash, 0);
  std::cout << "allocating closGeomGrids finished" << std::endl;

if( geom_hash == 1 )
{
  sim::Array<gitr_precision> hashX0(nHashes, 0.0), hashX1(nHashes, 0.0),
      hashY0(nHashes, 0.0), hashY1(nHashes, 0.0), hashZ0(nHashes, 0.0),
      hashZ1(nHashes, 0.0);
  if (world_rank == 0) {
    libconfig::Setting &geomHash = cfg.lookup("geometry_hash");
    if (nHashes > 1) {
      for (int i = 0; i < nHashes; i++) {
        hashX0[i] = geomHash["hashX0"][i];
        hashX1[i] = geomHash["hashX1"][i];
        hashZ0[i] = geomHash["hashZ0"][i];
        hashZ1[i] = geomHash["hashZ1"][i];
        if( use_3d_geom > 0 )
        {
          hashY0[i] = geomHash["hashY0"][i];
          hashY1[i] = geomHash["hashY1"][i];
        }
      }
    } else {
      getVariable(cfg, geomHashCfg + "hashX0", hashX0[0]);
      getVariable(cfg, geomHashCfg + "hashX1", hashX1[0]);
      getVariable(cfg, geomHashCfg + "hashZ0", hashZ0[0]);
      getVariable(cfg, geomHashCfg + "hashZ1", hashZ1[0]);
      if( use_3d_geom > 0 )
      {
      getVariable(cfg, geomHashCfg + "hashY0", hashY0[0]);
      getVariable(cfg, geomHashCfg + "hashY1", hashY1[0]);
      }
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&hashX0[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashX1[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashY0[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashY1[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashZ0[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashZ1[0], nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  int nHash = 0;
  int hashSum = 0;
  for (int i = 0; i < nR_closeGeomTotal; i++) {
    if (i == hashSum + nR_closeGeom[nHash]) {
      hashSum = hashSum + nR_closeGeom[nHash];
      nHash = nHash + 1;
    }
    closeGeomGridr[i] = (hashX1[nHash] - hashX0[nHash]) * (i - hashSum) /
                            (nR_closeGeom[nHash] - 1) +
                        hashX0[nHash];
  }
  nHash = 0;
  hashSum = 0;
  for (int j = 0; j < nY_closeGeomTotal; j++) {
    if (j == hashSum + nY_closeGeom[nHash]) {
      hashSum = hashSum + nY_closeGeom[nHash];
      nHash = nHash + 1;
    }
    closeGeomGridy[j] = (hashY1[nHash] - hashY0[nHash]) * (j - hashSum) /
                            (nY_closeGeom[nHash] - 1) +
                        hashY0[nHash];
  }
  nHash = 0;
  hashSum = 0;
  for (int k = 0; k < nZ_closeGeomTotal; k++) {
    if (k == hashSum + nZ_closeGeom[nHash]) {
      hashSum = hashSum + nZ_closeGeom[nHash];
      nHash = nHash + 1;
    }
    closeGeomGridz[k] = (hashZ1[nHash] - hashZ0[nHash]) * (k - hashSum) /
                            (nZ_closeGeom[nHash] - 1) +
                        hashZ0[nHash];
  }

  std::cout << "about to create iterator1 " << std::endl;
  thrust::counting_iterator<std::size_t> lines0(0);
  std::cout << "iterator2 " << std::endl;
  thrust::counting_iterator<std::size_t> lines1(nHashPointsTotal);
  int nHashMeshPointsPerProcess = ceil(nHashPointsTotal / world_size);
  std::cout << "nHashMeshPointsPerProcess " << nHashMeshPointsPerProcess
            << std::endl;
  std::vector<int> hashMeshIncrements(world_size);
  for (int j = 0; j < world_size - 1; j++) {
    hashMeshIncrements[j] = nHashMeshPointsPerProcess;
    std::cout << "hashMeshIncrements " << hashMeshIncrements[j] << std::endl;
  }
  hashMeshIncrements[world_size - 1] =
      nHashPointsTotal - (world_size - 1) * nHashMeshPointsPerProcess;
  std::cout << "minDist1 " << nGeomHash << std::endl;
  std::cout << "nHashPointsTotal " << nHashPointsTotal << std::endl;
  int Maxn_closeGeomElements = 0;
  for (int i = 0; i < nHashes; i++) {
    if (n_closeGeomElements[i] > Maxn_closeGeomElements) {
      Maxn_closeGeomElements = n_closeGeomElements[i];
    }
  }

  std::cout << "Maxn_closeGeomElements " << Maxn_closeGeomElements << std::endl;
  sim::Array<gitr_precision> minDist1(Maxn_closeGeomElements, 1e6);
  std::cout << "Generating geometry hash" << sizeof(int) << " bytes per int, "
            << nGeomHash << " for the entire hash " << std::endl;

#if USE_CUDA > 0
#endif
  std::cout << "starting geomhash" << std::endl;
  typedef std::chrono::high_resolution_clock Time0;
  typedef std::chrono::duration<gitr_precision> fsec0;
  auto start_clock0 = Time0::now();
  
  std::cout << "geo1 numbers " << nLines << " "
        << nHashes  << " " << nR_closeGeom[0] <<  " "
        << nY_closeGeom[0]  << " " << nZ_closeGeom[0] <<" "
        << n_closeGeomElements[0] << std::endl;
  
  hashGeom geo1(nLines, nHashes, boundaries.data(), closeGeomGridr.data(),
                closeGeomGridy.data(), closeGeomGridz.data(),
                n_closeGeomElements.data(), closeGeom.data(),
                nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(), use_3d_geom );
  std::cout << "nHashPoints start stop " << world_rank * nHashMeshPointsPerProcess << " "
        << world_rank * nHashMeshPointsPerProcess + hashMeshIncrements[world_rank] - 1<< std::endl;
  thrust::for_each(thrust::device,
                   lines0 + world_rank * nHashMeshPointsPerProcess,
                   lines0 + world_rank * nHashMeshPointsPerProcess +
                       hashMeshIncrements[world_rank],
                   geo1);
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "starting mpi close geom" << world_rank << std::endl;
  int hashPoint = 0;
  int closeGeomPoint = 0;
  int closeGeomPointTotal = 0;
  int closeGeomPointProcIndex = 0;
  std::vector<int> closeGeomPointIndex(world_size, 0);
  std::vector<int> closeGeomPointIncrements(world_size, 0);
  std::cout << "starting mpi close geom2" << world_rank << std::endl;
  if (world_rank == 0) {
    for (int ii = 0; ii < nHashes; ii++) {
      for (int i = 0; i < nR_closeGeom[ii]; i++) {
        /* Why is a lambda used here? Conditionals were previously preprocessor directives,
           and control flow cannot be replicated by just swapping them out with if statements */
        auto f_0 =
        [ & ]()-> void
        {
          for (int k = 0; k < nZ_closeGeom[ii]; k++) {
            if (hashPoint == hashMeshIncrements[closeGeomPointProcIndex]) {
              closeGeomPointIndex[closeGeomPointProcIndex + 1] =
                  closeGeomPointTotal;
              closeGeomPointIncrements[closeGeomPointProcIndex] =
                  closeGeomPoint;
              closeGeomPointProcIndex = closeGeomPointProcIndex + 1;
              hashPoint = 0;
              closeGeomPoint = 0;
            }
            hashPoint = hashPoint + 1;
            closeGeomPoint = closeGeomPoint + n_closeGeomElements[ii];
            closeGeomPointTotal = closeGeomPointTotal + n_closeGeomElements[ii];
          }
        };

        if(use_3d_geom)
        {
          for (int j = 0; j < nY_closeGeom[ii]; j++)
          {
            f_0();
          }
        }
      }
    }

    std::cout << "starting mpi close geom3" << world_rank << std::endl;
    closeGeomPointIncrements[world_size - 1] = closeGeomPoint;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&closeGeomPointIndex[0], world_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&closeGeomPointIncrements[0], world_size, MPI_INT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&closeGeomPointTotal, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "collect closegeom " << world_rank << std::endl;
  // Collect stuff
  for (int rr = 1; rr < world_size; rr++) {
    if (world_rank == rr) {
      std::cout << " node sending " << world_rank << " "
                << closeGeomPointIndex[rr] << " "
                << closeGeomPointIncrements[rr] << std::endl;
      MPI_Send(&closeGeom[closeGeomPointIndex[rr]],
               closeGeomPointIncrements[rr], MPI_INT, 0, 1, MPI_COMM_WORLD);
    } else if (world_rank == 0) {
      MPI_Recv(&closeGeom[closeGeomPointIndex[rr]],
               closeGeomPointIncrements[rr], MPI_INT, rr, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      std::cout << " node receiving " << world_rank << " "
                << closeGeomPointIndex[rr] << " "
                << closeGeomPointIncrements[rr] << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(closeGeom.data(), closeGeomPointTotal, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
  auto finish_clock0 = Time0::now();
  fsec0 fs0 = finish_clock0 - start_clock0;
  printf("Time taken          is %6.3f (secs) \n", fs0.count());
  if (world_rank == 0) {
    for (int i = 0; i < nHashes; i++) {
      std::cout << "opening file" << std::endl;
      netCDF::NcFile ncFile_hash("output/geomHash" + std::to_string(i) + ".nc",
                         netCDF::NcFile::replace);
      std::cout << "opened file" << std::endl;
      netCDF::NcDim hashNR = ncFile_hash.addDim("nR", nR_closeGeom[i]);
      netCDF::NcDim hashNY;
      if( use_3d_geom > 0 )
      {
        hashNY = ncFile_hash.addDim("nY", nY_closeGeom[i]);
      }
      netCDF::NcDim hashNZ = ncFile_hash.addDim("nZ", nZ_closeGeom[i]);
      netCDF::NcDim hashN = ncFile_hash.addDim("n", n_closeGeomElements[i]);
      vector<netCDF::NcDim> geomHashDim;
      geomHashDim.push_back(hashNR);
      std::cout << "created dims" << std::endl;
      if( use_3d_geom > 0 )
      {
        geomHashDim.push_back(hashNY);
      }
      geomHashDim.push_back(hashNZ);
      geomHashDim.push_back(hashN);
      netCDF::NcVar hash_gridR = ncFile_hash.addVar("gridR", netcdf_precision, hashNR);
      std::cout << "created dims2" << std::endl;
      netCDF::NcVar hash_gridY;
      if( use_3d_geom > 0 )
      {
        hash_gridY = ncFile_hash.addVar("gridY", netcdf_precision, hashNY);
      }
      netCDF::NcVar hash_gridZ = ncFile_hash.addVar("gridZ", netcdf_precision, hashNZ);
      netCDF::NcVar hash = ncFile_hash.addVar("hash", netCDF::ncInt, geomHashDim);
      std::cout << "created vars" << std::endl;
      int ncIndex = 0;
      if (i > 0)
        ncIndex = nR_closeGeom[i - 1];
      hash_gridR.putVar(&closeGeomGridr[ncIndex]);
      if( use_3d_geom > 0 )
      {
        if (i > 0)
          ncIndex = nY_closeGeom[i - 1];
        hash_gridY.putVar(&closeGeomGridy[ncIndex]);
      }

      if (i > 0)
        ncIndex = nZ_closeGeom[i - 1];
      hash_gridZ.putVar(&closeGeomGridz[ncIndex]);
      if (i > 0)
        ncIndex = nR_closeGeom[i - 1] * nY_closeGeom[i - 1] *
                  nZ_closeGeom[i - 1] * n_closeGeomElements[i - 1];
      hash.putVar(&closeGeom[ncIndex]);
      ncFile_hash.close();
    }
  }
      std::cout << "created vars2" << std::endl;
}

else if( geom_hash > 1 )
{
  if (world_rank == 0) {
    for (int i = 0; i < nHashes; i++) {
      int dataIndex = 0;
      if (i > 0)
        dataIndex = nR_closeGeom[0];
      getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg, "gridRString",
                     closeGeomGridr[dataIndex]);
      if (i > 0)
        dataIndex = nZ_closeGeom[0];
      getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg, "gridZString",
                     closeGeomGridz[dataIndex]);
      std::cout << "created vars3" << std::endl;
      if( use_3d_geom > 0 )
      {
        if (i > 0)
          dataIndex = nY_closeGeom[0];
        getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg, "gridYString",
                       closeGeomGridy[dataIndex]);
      }
      if (i > 0)
        dataIndex = nR_closeGeom[0] * nY_closeGeom[0] * nZ_closeGeom[0] *
                    n_closeGeomElements[0];
      getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg,
                     "closeGeomString", closeGeom[dataIndex]);
    }
  }
      std::cout << "created vars4" << std::endl;
#if USE_MPI > 0
  MPI_Bcast(closeGeomGridr.data(), nR_closeGeomTotal, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeomGridy.data(), nY_closeGeomTotal, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeomGridz.data(), nZ_closeGeomTotal, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeom.data(), nGeomHash, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

  int nR_closeGeom_sheath = 1;
  int nY_closeGeom_sheath = 1;
  int nZ_closeGeom_sheath = 1;
  int n_closeGeomElements_sheath = 1;
  int nGeomHash_sheath = 1;
  std::string geomHashSheathCfg = "geometry_sheath.";
if( geom_hash_sheath == 1 )
{
  if (world_rank == 0) {
    getVariable(cfg, geomHashSheathCfg + "nR_closeGeom", nR_closeGeom_sheath);
    getVariable(cfg, geomHashSheathCfg + "nZ_closeGeom", nZ_closeGeom_sheath);
    getVariable(cfg, geomHashSheathCfg + "n_closeGeomElements",
                n_closeGeomElements_sheath);
    nGeomHash_sheath =
        nR_closeGeom_sheath * nZ_closeGeom_sheath * n_closeGeomElements_sheath;
    if( use_3d_geom > 0 )
    {
      getVariable(cfg, geomHashSheathCfg + "nY_closeGeom", nY_closeGeom_sheath);
      nGeomHash_sheath = nY_closeGeom_sheath * nGeomHash_sheath;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

std::string hashFile_sheath;
if( geom_hash_sheath > 1 )
{
  if (world_rank == 0) {
    getVariable(cfg, geomHashSheathCfg + "fileString", hashFile_sheath);
    nR_closeGeom_sheath = getDimFromFile(cfg, input_path + hashFile_sheath,
                                         geomHashSheathCfg, "gridNrString");
    nZ_closeGeom_sheath = getDimFromFile(cfg, input_path + hashFile_sheath,
                                         geomHashSheathCfg, "gridNzString");
    n_closeGeomElements_sheath =
        getDimFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                       "nearestNelementsString");
    nGeomHash_sheath =
        nR_closeGeom_sheath * nZ_closeGeom_sheath * n_closeGeomElements_sheath;
    if( use_3d_geom > 0 )
    {
      nY_closeGeom_sheath = getDimFromFile(cfg, input_path + hashFile_sheath,
                                           geomHashSheathCfg, "gridNyString");
      nGeomHash_sheath = nY_closeGeom_sheath * nGeomHash_sheath;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

  sim::Array<gitr_precision> closeGeomGridr_sheath(nR_closeGeom_sheath),
      closeGeomGridy_sheath(nY_closeGeom_sheath),
      closeGeomGridz_sheath(nZ_closeGeom_sheath);
  sim::Array<int> closeGeom_sheath(nGeomHash_sheath);

  if( geom_hash_sheath == 1 )
  {
  gitr_precision hashX0_s, hashX1_s, hashY0_s, hashY1_s, hashZ0_s, hashZ1_s;
  if (world_rank == 0) {
    getVariable(cfg, geomHashSheathCfg + "hashX0", hashX0_s);
    getVariable(cfg, geomHashSheathCfg + "hashX1", hashX1_s);
    getVariable(cfg, geomHashSheathCfg + "hashZ0", hashZ0_s);
    getVariable(cfg, geomHashSheathCfg + "hashZ1", hashZ1_s);
    if( use_3d_geom > 0 )
    {
    getVariable(cfg, geomHashSheathCfg + "hashY0", hashY0_s);
    getVariable(cfg, geomHashSheathCfg + "hashY1", hashY1_s);
//#endif
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&hashX0_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashX1_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashY0_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashY1_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashZ0_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hashZ1_s, nHashes, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  for (int i = 0; i < nR_closeGeom_sheath; i++) {
    closeGeomGridr_sheath[i] =
        (hashX1_s - hashX0_s) * i / (nR_closeGeom_sheath - 1) + hashX0_s;
  }
  for (int j = 0; j < nY_closeGeom_sheath; j++) {
    closeGeomGridy_sheath[j] =
        (hashY1_s - hashY0_s) * j / (nY_closeGeom_sheath - 1) + hashY0_s;
  }
  for (int k = 0; k < nZ_closeGeom_sheath; k++) {
    closeGeomGridz_sheath[k] =
        (hashZ1_s - hashZ0_s) * k / (nZ_closeGeom_sheath - 1) + hashZ0_s;
  }

  thrust::counting_iterator<std::size_t> lines0_s(0);
  thrust::counting_iterator<std::size_t> lines1_s(nR_closeGeom_sheath *
                                                  nY_closeGeom_sheath);
  sim::Array<gitr_precision> minDist1_s(nGeomHash_sheath, 1e6);
  int nHashMeshPointsPerProcess_s =
      ceil(nR_closeGeom_sheath * nY_closeGeom_sheath * nZ_closeGeom_sheath /
           world_size);
  std::vector<int> hashMeshIncrements_s(world_size);
  for (int j = 0; j < world_size - 1; j++) {
    hashMeshIncrements_s[j] = nHashMeshPointsPerProcess_s;
  }
  hashMeshIncrements_s[world_size - 1] =
      nR_closeGeom_sheath * nY_closeGeom_sheath * nZ_closeGeom_sheath -
      (world_size - 1) * nHashMeshPointsPerProcess_s;
  typedef std::chrono::high_resolution_clock Time0_s;
  typedef std::chrono::duration<gitr_precision> fsec0_s;
  auto start_clock0_s = Time0_s::now();
  hashGeom_sheath geo_s(
      nLines, boundaries.data(), closeGeomGridr_sheath.data(),
      closeGeomGridy_sheath.data(), closeGeomGridz_sheath.data(),
      n_closeGeomElements_sheath, closeGeom_sheath.data(), nR_closeGeom_sheath,
      nY_closeGeom_sheath, nZ_closeGeom_sheath, use_3d_geom );
  thrust::for_each(thrust::device,
                   lines0_s + world_rank * nHashMeshPointsPerProcess_s,
                   lines0_s + world_rank * nHashMeshPointsPerProcess_s +
                       hashMeshIncrements_s[world_rank] - 1,
                   geo_s);
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
  MPI_Barrier(MPI_COMM_WORLD);
  // Collect stuff
  for (int rr = 1; rr < world_size; rr++) {
    if (world_rank == rr) {
      MPI_Send(&closeGeom_sheath[world_rank * nHashMeshPointsPerProcess_s *
                                 n_closeGeomElements_sheath],
               hashMeshIncrements_s[world_rank] * n_closeGeomElements_sheath,
               MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else if (world_rank == 0) {
      MPI_Recv(&closeGeom_sheath[rr * nHashMeshPointsPerProcess_s *
                                 n_closeGeomElements_sheath],
               hashMeshIncrements_s[rr] * n_closeGeomElements_sheath, MPI_INT,
               rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(closeGeom_sheath.data(),
            nR_closeGeom_sheath * nY_closeGeom_sheath * nZ_closeGeom_sheath *
                n_closeGeomElements_sheath,
            MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
  auto finish_clock0_s = Time0_s::now();
  fsec0_s fs0_s = finish_clock0_s - start_clock0_s;
  printf("Time taken          is %6.3f (secs) \n", fs0_s.count());
  if (world_rank == 0) {
    netCDF::NcFile ncFile_hash_sheath("output/geomHash_sheath.nc", netCDF::NcFile::replace);
    netCDF::NcDim hashNR_sheath = ncFile_hash_sheath.addDim("nR", nR_closeGeom_sheath);
    netCDF::NcDim hashNY_sheath = ncFile_hash_sheath.addDim("nY", nY_closeGeom_sheath);
    netCDF::NcDim hashNZ_sheath = ncFile_hash_sheath.addDim("nZ", nZ_closeGeom_sheath);
    netCDF::NcDim hashN_sheath =
        ncFile_hash_sheath.addDim("n", n_closeGeomElements_sheath);
    vector<netCDF::NcDim> geomHashDim_sheath;
    geomHashDim_sheath.push_back(hashNR_sheath);
    geomHashDim_sheath.push_back(hashNY_sheath);
    geomHashDim_sheath.push_back(hashNZ_sheath);
    geomHashDim_sheath.push_back(hashN_sheath);
    netCDF::NcVar hash_gridR_sheath =
        ncFile_hash_sheath.addVar("gridR", netcdf_precision, hashNR_sheath);
    netCDF::NcVar hash_gridY_sheath =
        ncFile_hash_sheath.addVar("gridY", netcdf_precision, hashNY_sheath);
    netCDF::NcVar hash_gridZ_sheath =
        ncFile_hash_sheath.addVar("gridZ", netcdf_precision, hashNZ_sheath);
    netCDF::NcVar hash_sheath =
        ncFile_hash_sheath.addVar("hash", netCDF::ncInt, geomHashDim_sheath);
    hash_gridR_sheath.putVar(&closeGeomGridr_sheath[0]);
    hash_gridY_sheath.putVar(&closeGeomGridy_sheath[0]);
    hash_gridZ_sheath.putVar(&closeGeomGridz_sheath[0]);
    hash_sheath.putVar(&closeGeom_sheath[0]);
    ncFile_hash_sheath.close();
  }
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
}
else if( geom_hash_sheath > 1 )
{
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridRString", closeGeomGridr_sheath[0]);
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridZString", closeGeomGridz_sheath[0]);
    if( use_3d_geom > 0 )
    {
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridYString", closeGeomGridy_sheath[0]);
//#endif
    }
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "closeGeomString", closeGeom_sheath[0]);
#if USE_MPI > 0
  }
  MPI_Bcast(closeGeomGridr_sheath.data(), nR_closeGeom_sheath, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeomGridy_sheath.data(), nY_closeGeom_sheath, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeomGridz_sheath.data(), nZ_closeGeom_sheath, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(closeGeom_sheath.data(), nGeomHash_sheath, MPI_INT, 0,
            MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

  int nR_Lc = 1;
  int nY_Lc = 1;
  int nZ_Lc = 1;
  int nTracers = 1;
  gitr_precision r0_Lc, r1_Lc, y0_Lc, y1_Lc, z0_Lc, z1_Lc, dr;
  int nTraceSteps;
  std::string connLengthCfg = "connectionLength.";
  std::string lcFile;
    NcFile ncFileLC("output/LcS.nc", NcFile::replace);
    NcDim nc_nTracers = ncFileLC.addDim("nTracers", nTracers);

if (world_rank == 0) {
if( GENERATE_LC > 0 )
{
    getVariable(cfg, connLengthCfg + "fileString", lcFile);
    getVariable(cfg, connLengthCfg + "nX", nR_Lc);
    getVariable(cfg, connLengthCfg + "nY", nY_Lc);
    getVariable(cfg, connLengthCfg + "nZ", nZ_Lc);
    getVariable(cfg, connLengthCfg + "netx0", r0_Lc);
    getVariable(cfg, connLengthCfg + "netx1", r1_Lc);
    getVariable(cfg, connLengthCfg + "nety0", y0_Lc);
    getVariable(cfg, connLengthCfg + "nety1", y1_Lc);
    getVariable(cfg, connLengthCfg + "netz0", z0_Lc);
    getVariable(cfg, connLengthCfg + "netz1", z1_Lc);
    getVariable(cfg, connLengthCfg + "nTraceSteps", nTraceSteps);
    getVariable(cfg, connLengthCfg + "dr", dr);
}
else
{
  if( LC_INTERP > 1 )
  {
      nR_Lc =
          getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNrString");
      nY_Lc =
          getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNyString");
      nZ_Lc =
          getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNzString");
  }
}
}
#if USE_MPI > 0
  MPI_Bcast(&nR_Lc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Lc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Lc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r0_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r1_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y0_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y1_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&z0_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&z1_Lc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nTraceSteps, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dr, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

if( use_3d_geom > 0 )
{
  nTracers = nR_Lc * nY_Lc * nZ_Lc;
}
else
{
  nTracers = nR_Lc * nZ_Lc;
}

  sim::Array<gitr_precision> Lc(nTracers), s(nTracers);
  sim::Array<gitr_precision> gridRLc(nR_Lc), gridYLc(nY_Lc), gridZLc(nZ_Lc);
  sim::Array<int> noIntersectionNodes(nTracers);

if( GENERATE_LC > 0 )
{
  gitr_precision lcBuffer = 0.0;
    if( use_3d_geom > 0 )
    {
  gitr_precision dy_Lc = (y1_Lc - y0_Lc) / (nY_Lc - 1);
  for (int j = 0; j < nY_Lc; j++) {
    gridYLc[j] = y0_Lc + j * dy_Lc;
  }
//#endif
    }
  gitr_precision dr_Lc = (r1_Lc - r0_Lc) / (nR_Lc - 1);
  for (int i = 0; i < nR_Lc; i++) {
    gridRLc[i] = r0_Lc + i * dr_Lc;
  }

  gitr_precision dz_Lc = (z1_Lc - z0_Lc) / (nZ_Lc - 1);
  for (int j = 0; j < nZ_Lc; j++) {
    gridZLc[j] = z0_Lc + j * dz_Lc;
  }
  if (nTracers % world_size != 0) {
    std::cout << "nTracers for Lc not divisible by num MPI threads"
              << std::endl;
    exit(0);
  }
  // for(size_t i = 0; i < world_size; i++) {
  std::pair<size_t, size_t> pair{world_rank * nTracers / world_size,
                                 (world_rank + 1) * nTracers / world_size};
  std::cout << "Group " << (world_rank + 1) << ": [" << pair.first << ","
            << pair.second << ") - " << (pair.second - pair.first)
            << " elements." << std::endl;
  //}
  std::cout << "Creating tracer particles" << std::endl;
  thrust::counting_iterator<std::size_t> lcBegin(pair.first);
  thrust::counting_iterator<std::size_t> lcEnd(pair.second - 1);
  auto forwardTracerParticles = new Particles(nTracers, 1, cfg, gitr_flags );
  auto backwardTracerParticles = new Particles(nTracers, 1, cfg, gitr_flags );
  int addIndex = 0;
  std::cout << "Initializing tracer particles" << std::endl;

  for (int i = 0; i < nR_Lc; i++) {
    for (int j = 0; j < nY_Lc; j++) {
      for (int k = 0; k < nZ_Lc; k++) {
    if( use_3d_geom > 0 )
    {
        addIndex = i + j * nR_Lc + k * nR_Lc * nY_Lc;
    }
    else
    {
//#else
        addIndex = i + k * nR_Lc;
//#endif
    }
        forwardTracerParticles->setParticle(addIndex, gridRLc[i], gridYLc[j],
                                            gridZLc[k], 0.0, 0.0, 0.0, 0, 0.0,
                                            0.0, 1);
        backwardTracerParticles->setParticle(addIndex, gridRLc[i], gridYLc[j],
                                             gridZLc[k], 0.0, 0.0, 0.0, 0, 0.0,
                                             0.0,1);
      }
    }
  }

  // dummy surfaces for Lcs calculation (geometry_check)
  auto dummy_surfaces = new Surfaces(1, 1, 1, 1);
  dummy_surfaces->setSurface(1, 1, 1, 1, 1, 1);

  typedef std::chrono::high_resolution_clock Time_trace;
  typedef std::chrono::duration<gitr_precision> fsec_trace;
  auto start_clock_trace = Time_trace::now();
  std::cout << "Starting trace loop" << std::endl;
  std::cout << "nTraceSteps" << nTraceSteps << " dr " << dr << std::endl;
  for (int ii = 0; ii < nTraceSteps; ii++) {
#if USE_CUDA
    cudaDeviceSynchronize();
#endif
    thrust::for_each(thrust::device, lcBegin, lcEnd,
                     field_line_trace(1.0, forwardTracerParticles, dr,
                                      boundaries.data(), nLines, nR_Lc, nZ_Lc,
                                      gridRLc.data(), gridZLc.data(), Lc.data(),
                                      nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                                      &bfieldGridz.front(), &br.front(),
                                      &bz.front(), &by.front(), cylsymm ));

    thrust::for_each(thrust::device, lcBegin, lcEnd,
                     field_line_trace(-1.0, backwardTracerParticles, dr,
                                      boundaries.data(), nLines, nR_Lc, nZ_Lc,
                                      gridRLc.data(), gridZLc.data(), Lc.data(),
                                      nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                                      &bfieldGridz.front(), &br.front(),
                                      &bz.front(), &by.front(), cylsymm ));

    thrust::for_each(
        thrust::device, lcBegin, lcEnd,
        geometry_check(forwardTracerParticles, nLines, &boundaries[0],
                       dummy_surfaces, dr, ii, nR_closeGeom.data(),
                       nY_closeGeom.data(), nZ_closeGeom.data(),
                       n_closeGeomElements.data(), &closeGeomGridr.front(),
                       &closeGeomGridy.front(), &closeGeomGridz.front(),
                       &closeGeom.front(), 0, 0.0, 0.0, 0, 0.0, 0.0, flux_ea, surface_model,
                       geom_hash,
                       use_3d_geom,
                       cylsymm ) );

    thrust::for_each(
        thrust::device, lcBegin, lcEnd,
        geometry_check(backwardTracerParticles, nLines, &boundaries[0],
                       dummy_surfaces, dr, ii, nR_closeGeom.data(),
                       nY_closeGeom.data(), nZ_closeGeom.data(),
                       n_closeGeomElements.data(), &closeGeomGridr.front(),
                       &closeGeomGridy.front(), &closeGeomGridz.front(),
                       &closeGeom.front(), 0, 0.0, 0.0, 0, 0.0, 0.0, flux_ea, surface_model,
                       geom_hash,
                       use_3d_geom,
                       cylsymm ) );
  }
  auto finish_clock_trace = Time_trace::now();
  fsec_trace fstrace = finish_clock_trace - start_clock_trace;
  printf("Time taken          is %6.3f (secs) \n", fstrace.count());
  printf("Time taken per step is %6.3f (secs) \n",
         fstrace.count() / (gitr_precision)nTraceSteps);
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
  sim::Array<gitr_precision> forwardHitWall(nTracers, 0.0),
      backwardHitWall(nTracers, 0.0), forwardTracerX(nTracers, 0.0),
      backwardTracerX(nTracers, 0.0);
  sim::Array<gitr_precision> forwardTracerY(nTracers, 0.0),
      backwardTracerY(nTracers, 0.0);
  sim::Array<gitr_precision> forwardTracerZ(nTracers, 0.0),
      backwardTracerZ(nTracers, 0.0);
  sim::Array<gitr_precision> forwardDistanceTraveled(nTracers, 0.0),
      backwardDistanceTraveled(nTracers, 0.0);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(
      &forwardTracerParticles->hitWall[world_rank * nTracers / world_size],
      nTracers / world_size, MPI_FLOAT, &forwardHitWall[0],
      nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(
      &backwardTracerParticles->hitWall[world_rank * nTracers / world_size],
      nTracers / world_size, MPI_FLOAT, &backwardHitWall[0],
      nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&forwardTracerParticles
                  ->distanceTraveled[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &forwardDistanceTraveled[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&backwardTracerParticles
                  ->distanceTraveled[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &backwardDistanceTraveled[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&forwardTracerParticles->x[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &forwardTracerX[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&backwardTracerParticles->x[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &backwardTracerX[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&forwardTracerParticles->y[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &forwardTracerY[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&backwardTracerParticles->y[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &backwardTracerY[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&forwardTracerParticles->z[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &forwardTracerZ[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&backwardTracerParticles->z[world_rank * nTracers / world_size],
             nTracers / world_size, MPI_FLOAT, &backwardTracerZ[0],
             nTracers / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i = 0; i < nTracers; i++) {
    forwardTracerParticles->hitWall[i] = forwardHitWall[i];
    forwardTracerParticles->distanceTraveled[i] = forwardDistanceTraveled[i];
    backwardTracerParticles->hitWall[i] = backwardHitWall[i];
    backwardTracerParticles->distanceTraveled[i] = backwardDistanceTraveled[i];
  }
  if (world_rank == 0) {
#endif
    addIndex = 0;
    gitr_precision forwardDist = 0.0;
    gitr_precision backwardDist = 0.0;
    for (int i = 0; i < nR_Lc; i++) {
      for (int j = 0; j < nY_Lc; j++) {
        for (int k = 0; k < nZ_Lc; k++) {
    if( use_3d_geom > 0 )
    {
          addIndex = i + j * nR_Lc + k * nR_Lc * nY_Lc;
    }
    else
    {
        addIndex = i + k * nR_Lc;
    }
          if (forwardTracerParticles->hitWall[addIndex] > 0) {
            forwardDist = forwardTracerParticles->distanceTraveled[addIndex];
          } else {
            forwardDist = 0.0;
          }

          if (backwardTracerParticles->hitWall[addIndex] > 0) {
            backwardDist = backwardTracerParticles->distanceTraveled[addIndex];
          } else
            backwardDist = 0.0;

          Lc[addIndex] = forwardDist + backwardDist;
          if (Lc[addIndex] > 0.0) {
            Lc[addIndex] = Lc[addIndex] + lcBuffer;
          } else
            Lc[addIndex] = 1.0e4;
          if (forwardTracerParticles->hitWall[addIndex] > 0 &&
              backwardTracerParticles->hitWall[addIndex] > 0) {
            s[addIndex] =
                min(forwardTracerParticles->distanceTraveled[addIndex],
                    backwardTracerParticles->distanceTraveled[addIndex]);
            if (forwardTracerParticles->distanceTraveled[addIndex] >
                0.5 * Lc[addIndex]) {
              s[addIndex] = -s[addIndex];
            }
          } else if (forwardTracerParticles->hitWall[addIndex] > 0) {
            s[addIndex] = forwardTracerParticles->distanceTraveled[addIndex];
          } else if (backwardTracerParticles->hitWall[addIndex] > 0) {
            s[addIndex] = -backwardTracerParticles->distanceTraveled[addIndex];
          } else
            s[addIndex] = 1.0e4;
          if (forwardTracerParticles->hitWall[addIndex] +
                  backwardTracerParticles->hitWall[addIndex] <
              4.0) {
            noIntersectionNodes[addIndex] = 1;
          }
        }
      }
    }

    vector<NcDim> dims_lc;
    NcDim nc_nRLc = ncFileLC.addDim("nR", nR_Lc);
    dims_lc.push_back(nc_nRLc);

    NcDim nc_nYLc;
    if( use_3d_geom )
    {
    nc_nYLc = ncFileLC.addDim("nY", nY_Lc);
    dims_lc.push_back(nc_nYLc);
    }

    NcDim nc_nZLc = ncFileLC.addDim("nZ", nZ_Lc);
    dims_lc.push_back(nc_nZLc);

    NcVar nc_Lc = ncFileLC.addVar("Lc", netcdf_precision, dims_lc);
    NcVar nc_s = ncFileLC.addVar("s", netcdf_precision, dims_lc);
    NcVar nc_ftx = ncFileLC.addVar("fx", netcdf_precision, dims_lc);
    NcVar nc_fty = ncFileLC.addVar("fy", netcdf_precision, dims_lc);
    NcVar nc_ftz = ncFileLC.addVar("fz", netcdf_precision, dims_lc);
    NcVar nc_btx = ncFileLC.addVar("bx", netcdf_precision, dims_lc);
    NcVar nc_bty = ncFileLC.addVar("by", netcdf_precision, dims_lc);
    NcVar nc_btz = ncFileLC.addVar("bz", netcdf_precision, dims_lc);
    NcVar nc_nI = ncFileLC.addVar("noIntersection", netcdf_precision, dims_lc);
    NcVar nc_gridRLc = ncFileLC.addVar("gridR", netcdf_precision, nc_nRLc);
    NcVar nc_gridYLc;
    if( use_3d_geom )
    {
    nc_gridYLc = ncFileLC.addVar("gridY", netcdf_precision, nc_nYLc);
    }
    NcVar nc_gridZLc = ncFileLC.addVar("gridZ", netcdf_precision, nc_nZLc);

    if( use_3d_geom )
    {
    nc_gridYLc.putVar(&gridYLc[0]);
//#endif
    }
    nc_gridZLc.putVar(&gridZLc[0]);
    ncFileLC.close();
#if USE_MPI > 0
  }
  MPI_Bcast(Lc.data(), nTracers, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(s.data(), nTracers, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(noIntersectionNodes.data(), nTracers, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
  //}
}

if( LC_INTERP > 1 )
{
  std::cout << "Importing pre-existing connection length file" << std::endl;
  getVariable(cfg, connLengthCfg + "fileString", lcFile);
  getVarFromFile(cfg, input_path + lcFile, connLengthCfg, "gridRString",
                 gridRLc);
  getVarFromFile(cfg, input_path + lcFile, connLengthCfg, "gridYString",
                 gridYLc);
  getVarFromFile(cfg, input_path + lcFile, connLengthCfg, "gridZString",
                 gridZLc);
  getVarFromFile(cfg, input_path + lcFile, connLengthCfg, "LcString", Lc);
  getVarFromFile(cfg, input_path + lcFile, connLengthCfg, "SString", s);
}


  // Background Plasma Temperature Initialization
  int nR_Temp = 1;
  int nY_Temp = 1;
  int nZ_Temp = 1;
  int n_Temp = 1;
  std::string tempCfg = "backgroundPlasmaProfiles.Temperature.";
  std::string tempFile;
if(temp_interp > 0 )
{
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, tempCfg + "fileString", tempFile);
    nR_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNrString");
        if(temp_interp > 1 )
        {
    nZ_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNzString");
        }
        if( temp_interp > 2 )
        {
    nY_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNyString");
        }
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


  sim::Array<gitr_precision> TempGridr(nR_Temp), TempGridz(nZ_Temp), TempGridy(nY_Temp);
  n_Temp = nR_Temp * nY_Temp * nZ_Temp;
  sim::Array<gitr_precision> ti(n_Temp), te(n_Temp);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    if( temp_interp == 0 )
    {
    getVariable(cfg, tempCfg + "ti", ti[0]);
    getVariable(cfg, tempCfg + "te", te[0]);
    }
    else
    {
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridRString",
                 TempGridr[0]);
  if( temp_interp > 1 )
  {
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridZString",
                 TempGridz[0]);
  }
  if( temp_interp > 2 )
  {
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridYString",
                 TempGridy[0]);
  }
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "IonTempString", ti[0]);
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "ElectronTempString",
                 te[0]);
    }
#if USE_MPI > 0
  }

  MPI_Bcast(TempGridr.data(), nR_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(TempGridy.data(), nY_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(TempGridz.data(), nZ_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ti.data(), n_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(te.data(), n_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  gitr_precision testVec = 0.0;
  testVec = interp2dCombined(0.0, 0.1, 0.0, nR_Temp, nZ_Temp, TempGridr.data(),
                             TempGridz.data(), ti.data(), cylsymm );
  std::cout << "Finished Temperature import " << testVec << std::endl;

  // Background Plasma Density Initialization
  int nR_Dens = 1;
  int nY_Dens = 1;
  int nZ_Dens = 1;
  int n_Dens = 1;
  std::string densCfg = "backgroundPlasmaProfiles.Density.";
  std::string densFile;
  if( density_interp > 0 )
  {
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, densCfg + "fileString", densFile);
    nR_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNrString");
        if( density_interp > 1 )
        {
    nZ_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNzString");
        }
        if(density_interp > 2 )
        {
    nY_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNyString");
        }
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  sim::Array<gitr_precision> DensGridr(nR_Dens), DensGridz(nZ_Dens), DensGridy(nY_Dens);
  n_Dens = nR_Dens * nY_Dens * nZ_Dens;
  sim::Array<gitr_precision> ni(n_Dens), ne(n_Dens);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
if( density_interp == 0 )
{
    getVariable(cfg, densCfg + "ni", ni[0]);
    getVariable(cfg, densCfg + "ne", ne[0]);
    } else {
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridRString",
                 DensGridr[0]);
  if( density_interp > 1 )
  {
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridZString",
                 DensGridz[0]);
  }
  if( density_interp > 2 )
  {
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridYString",
                 DensGridy[0]);
  }
  getVarFromFile(cfg, input_path + densFile, densCfg, "IonDensityString",
                 ni[0]);
  getVarFromFile(cfg, input_path + densFile, densCfg, "ElectronDensityString",
                 ne[0]);
}
    std::cout << "Finished density import "
              << interp2dCombined(5.5, 0.0, -4.4, nR_Dens, nZ_Dens,
                                  &DensGridr.front(), &DensGridz.front(),
                                  &ne.front(), cylsymm )
              << " "
              << interp2dCombined(0.0, 0.1, 0.0, nR_Dens, nZ_Dens,
                                  &DensGridr.front(), &DensGridz.front(),
                                  &ne.front(), cylsymm )
              << std::endl;

#if USE_MPI > 0
  }
  MPI_Bcast(DensGridr.data(), nR_Dens, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(DensGridy.data(), nY_Dens, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(DensGridz.data(), nZ_Dens, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ni.data(), n_Dens, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ne.data(), n_Dens, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  // Background Plasma flow velocity initialization
  std::cout << "Starting flow import " << endl;
  int nR_flowV = 1;
  int nY_flowV = 1;
  int nZ_flowV = 1;
  int n_flowV = 1;
  std::string flowVCfg = "backgroundPlasmaProfiles.FlowVelocity.";
  if( flowv_interp == 1 )
  {
  nR_flowV = nR_Lc;
  nY_flowV = nY_Lc;
  nZ_flowV = nZ_Lc;
  }
  std::string flowVFile;
  if( flowv_interp > 1 )
  {
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, flowVCfg + "fileString", flowVFile);
    nR_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNrString");
    nZ_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNzString");
    if( flowv_interp > 2 )
    {
    nY_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNyString");
    }
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
  sim::Array<gitr_precision> flowVGridr(nR_flowV), flowVGridy(nY_flowV),
      flowVGridz(nZ_flowV);
  n_flowV = nR_flowV * nY_flowV * nZ_flowV;
  sim::Array<gitr_precision> flowVr(n_flowV), flowVz(n_flowV), flowVt(n_flowV);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
  if( flowv_interp == 0 )
  {
    getVariable(cfg, flowVCfg + "flowVr", flowVr[0]);
    getVariable(cfg, flowVCfg + "flowVy", flowVt[0]);
    getVariable(cfg, flowVCfg + "flowVz", flowVz[0]);
  }
  else
  {
    if( flowv_interp > 1 )
    {
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "gridRString",
                 flowVGridr[0]);
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "gridZString",
                 flowVGridz[0]);
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "flowVrString",
                 flowVr[0]);
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "flowVtString",
                 flowVt[0]);
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "flowVzString",
                 flowVz[0]);
    }
    if( flowv_interp > 2 )
    {
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "gridYString",
                 flowVGridy[0]);
    }
  }
#if USE_MPI > 0
  }
  MPI_Bcast(flowVGridr.data(), nR_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(flowVGridy.data(), nY_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(flowVGridz.data(), nZ_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(flowVr.data(), n_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(flowVt.data(), n_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(flowVz.data(), n_flowV, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  gitr_precision surroundingMinimumR = 0.0;
  gitr_precision surroundingMinimumY = 0.0;
  gitr_precision surroundingMinimumZ = 0.0;
  int iterIndex = 0;
  int nFlowVs;

if( flowv_interp == 1 )
{
  for (int i = 0; i < nR_flowV; i++) {
    std::cout << " !!! gridRLc " << gridRLc[i] << std::endl;
  }
  std::cout << " !!! gridZLc " << gridZLc[0] << std::endl;
  for (int i = 0; i < nR_flowV; i++) {
    flowVGridr[i] = gridRLc[i];
  }
  for (int i = 0; i < nZ_flowV; i++) {
    flowVGridz[i] = gridZLc[i];
  }
  std::cout << " !!! flowvgridr0 " << flowVGridr[0] << std::endl;
  nFlowVs = nR_Lc * nZ_Lc;
  if( LC_INTERP == 3 )
  {
  for (int i = 0; i < nY_flowV; i++)
    flowVGridy[i] = gridYLc[i];
  nFlowVs = nR_Lc * nY_Lc * nZ_Lc;
  }
  gitr_precision thisY = 0.0;
  gitr_precision cs0 = 0.0;
  gitr_precision teLocal = 0.0;
  gitr_precision tiLocal = 0.0;
  gitr_precision BLocal[3] = {0.0, 0.0, 0.0};
  gitr_precision Bnorm[3] = {0.0, 0.0, 0.0};
  gitr_precision Bmag = 0.0;
  int index = 0;
  gitr_precision cs = 0.0;
  gitr_precision absS = 0.0;
  std::cout << "Beginning analytic flowV calculation " << std::endl;
  for (int i = 0; i < nR_Lc; i++) {
#if LC_INTERP == 3
    for (int k = 0; k < nY_Lc; k++) {
      thisY = flowVGridy[k];
#endif

      for (int j = 0; j < nZ_Lc; j++) {
        // std::cout << "debug here 1 " << i << " " << j << std::endl;
        // std::cout << "debug here 2 " << flowVGridr[i] << " " << thisY << " "
        //    << flowVGridz[j] << " " << nR_Temp << " "<<nZ_Temp << std::endl;
        teLocal = interp2dCombined(flowVGridr[i], thisY, flowVGridz[j], nR_Temp,
                                   nZ_Temp, &TempGridr.front(),
                                   &TempGridz.front(), &te.front(), cylsymm );
        tiLocal = interp2dCombined(flowVGridr[i], thisY, flowVGridz[j], nR_Temp,
                                   nZ_Temp, &TempGridr.front(),
                                   &TempGridz.front(), &ti.front(), cylsymm );
        cs0 =
            std::sqrt((teLocal + tiLocal) * 1.602e-19 / (background_amu * 1.66e-27));
        interp2dVector(&BLocal[0], flowVGridr[i], thisY, flowVGridz[j],
                       nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                       bfieldGridz.data(), br.data(), bz.data(), by.data(), cylsymm );
        Bmag = std::sqrt(BLocal[0] * BLocal[0] + BLocal[1] * BLocal[1] +
                    BLocal[2] * BLocal[2]);
        Bnorm[0] = BLocal[0] / Bmag;
        Bnorm[1] = BLocal[1] / Bmag;
        Bnorm[2] = BLocal[2] / Bmag;

#if LC_INTERP == 3
        index = i + k * nR_Lc + j * nR_Lc * nY_Lc;
        // std::cout << "flowv calc index " << index << std::endl;
#else
      index = i + j * nR_Lc;
#endif
        absS = std::abs(s[index]);
        cs = cs0 * (0.5 * Lc[index] / absS -
                    std::sqrt(0.25 * Lc[index] * Lc[index] / absS / absS - 1.0));
        if (std::isnan(cs))
          cs = 0.0;
        flowVr[index] = std::copysign(1.0,s[index]) * Bnorm[0] * cs;
        flowVt[index] = std::copysign(1.0,s[index]) * Bnorm[1] * cs;
        flowVz[index] = std::copysign(1.0,s[index]) * Bnorm[2] * cs;
#if LC_INTERP == 3
      }
#endif
    }
  }
  std::cout << "Done with initial calculation, beginning sorting" << std::endl;
  sim::Array<gitr_precision> flowVrSub(nFlowVs), flowVzSub(nFlowVs), flowVySub(nFlowVs);
  sim::Array<int> noIntersectionNearestMax(nFlowVs);
  for (int i = 0; i < nR_Lc; i++) {
    std::cout << "i of " << i << " " << nR_Lc << std::endl;
    for (int j = 0; j < nY_Lc; j++) {
      for (int k = 0; k < nZ_Lc; k++) {
        index = i + j * nR_Lc + k * nR_Lc * nY_Lc;
        if (noIntersectionNodes[index] == 1) {
          surroundingMinimumR = 0.0;
          surroundingMinimumY = 0.0;
          surroundingMinimumZ = 0.0;
          for (int ii = i - 1; ii < i + 2; ii++) {
            for (int jj = j - 1; jj < j + 2; jj++) {
              for (int kk = k - 1; kk < k + 2; kk++) {
                iterIndex = ii + jj * nR_Lc + kk * nR_Lc * nY_Lc;
                if (iterIndex > 0 && iterIndex < nFlowVs) {
                  if (noIntersectionNodes[iterIndex] == 0) {
                    if (std::abs(flowVr[iterIndex]) > std::abs(surroundingMinimumR)) {
                      surroundingMinimumR = flowVr[iterIndex];
                    }
                    if (std::abs(flowVt[iterIndex]) > std::abs(surroundingMinimumY)) {
                      surroundingMinimumY = flowVt[iterIndex];
                    }
                    if (std::abs(flowVz[iterIndex]) > std::abs(surroundingMinimumZ)) {
                      surroundingMinimumZ = flowVz[iterIndex];
                    }
                  }
                }
              }
            }
          }
          flowVrSub[index] = surroundingMinimumR;
          flowVySub[index] = surroundingMinimumY;
          flowVzSub[index] = surroundingMinimumZ;
        }
      }
    }
  }
  for (int i = 0; i < nFlowVs; i++) {
    if (i == 282839) {
      std::cout << " noIntersectionNodes " << noIntersectionNodes[i]
                << std::endl;
    }
    if (noIntersectionNodes[i] == 1) {
      flowVr[i] = flowVrSub[i];
      flowVt[i] = flowVySub[i];
      flowVz[i] = flowVzSub[i];
    }
  }
  NcFile ncFileFlow("flowV.nc", NcFile::replace);
  NcDim nFlowV = ncFileFlow.addDim("n_flowV", n_flowV);
  NcDim nc_nRflow = ncFileFlow.addDim("nR", nR_flowV);
  NcDim nc_nYflow = ncFileFlow.addDim("nY", nY_flowV);
  NcDim nc_nZflow = ncFileFlow.addDim("nZ", nZ_flowV);
  vector<NcDim> dimsFlowV;
  dimsFlowV.push_back(nc_nZflow);
  dimsFlowV.push_back(nc_nYflow);
  dimsFlowV.push_back(nc_nRflow);
  NcVar nc_flowVr = ncFileFlow.addVar("flowVr", netcdf_precision, dimsFlowV);
  NcVar nc_flowVt = ncFileFlow.addVar("flowVt", netcdf_precision, dimsFlowV);
  NcVar nc_flowVz = ncFileFlow.addVar("flowVz", netcdf_precision, dimsFlowV);
  nc_flowVr.putVar(&flowVr[0]);
  nc_flowVt.putVar(&flowVt[0]);
  nc_flowVz.putVar(&flowVz[0]);
  ncFileFlow.close();
  std::string outnameFlowVr = "flowVr.m";
  std::string outnameFlowVz = "flowVz.m";
  std::string outnameFlowVt = "flowVt.m";
#if LC_INTERP == 3
  OUTPUT3d(profiles_folder, outnameFlowVr, nR_flowV, nY_flowV, nZ_flowV,
           &flowVr.front());
  OUTPUT3d(profiles_folder, outnameFlowVz, nR_flowV, nY_flowV, nZ_flowV,
           &flowVz.front());
  OUTPUT3d(profiles_folder, outnameFlowVt, nR_flowV, nY_flowV, nZ_flowV,
           &flowVt.front());
#else
  OUTPUT2d(profiles_folder, outnameFlowVr, nR_flowV, nZ_flowV, &flowVr.front());
  OUTPUT2d(profiles_folder, outnameFlowVz, nR_flowV, nZ_flowV, &flowVz.front());
  OUTPUT2d(profiles_folder, outnameFlowVt, nR_flowV, nZ_flowV, &flowVt.front());
#endif
}

  // Background plasma temperature gradient field intitialization
  int nR_gradT = 1;
  int nY_gradT = 1;
  int nZ_gradT = 1;
  int n_gradT = 1;
  std::string gradTCfg = "backgroundPlasmaProfiles.gradT.";
  std::string gradTFile;
  if (world_rank == 0) {
    if( gradt_interp > 0 )
    {
    getVariable(cfg, gradTCfg + "fileString", gradTFile);
    nR_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNrString");
    }

    if( gradt_interp > 1 )
    {
    nZ_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNzString");
    }

    if( gradt_interp > 2 )
    {
    nY_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNyString");
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  n_gradT = nR_gradT * nY_gradT * nZ_gradT;
  sim::Array<gitr_precision> gradTGridr(nR_gradT), gradTGridy(nY_gradT),
      gradTGridz(nZ_gradT);
  sim::Array<gitr_precision> gradTeR(n_gradT), gradTeZ(n_gradT), gradTeY(n_gradT),
      gradTiR(n_gradT), gradTiZ(n_gradT), gradTiY(n_gradT);

  if (world_rank == 0) {

    if( gradt_interp == 0 )
    {
    getVariable(cfg, gradTCfg + "gradTeR", gradTeR[0]);
    getVariable(cfg, gradTCfg + "gradTeY", gradTeY[0]);
    getVariable(cfg, gradTCfg + "gradTeZ", gradTeZ[0]);
    getVariable(cfg, gradTCfg + "gradTiR", gradTiR[0]);
    getVariable(cfg, gradTCfg + "gradTiY", gradTiY[0]);
    getVariable(cfg, gradTCfg + "gradTiZ", gradTiZ[0]);
    }
    else
    {
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridRString",
                   gradTGridr[0]);
    if( gradt_interp > 1 )
    {
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridZString",
                   gradTGridz[0]);
    }

    if( gradt_interp > 2 )
    {
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridYString",
                   gradTGridy[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTeYString",
                   gradTeY[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTiYString",
                   gradTiY[0]);
    }

    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTiRString",
                   gradTiR[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTiZString",
                   gradTiZ[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTeRString",
                   gradTeR[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTeZString",
                   gradTeZ[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTeYString",
                   gradTeY[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTiYString",
                   gradTiY[0]);
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&gradTGridr[0], nR_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTGridy[0], nY_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTGridz[0], nZ_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTeR[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTiR[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTeY[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTiY[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTeZ[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gradTiZ[0], n_gradT, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  gitr_precision gradTi[3] = {0.0};
  interp2dVector(&gradTi[0], 1.45, 0.0, -1.2, nR_gradT, nZ_gradT,
                 gradTGridr.data(), gradTGridz.data(), gradTiR.data(),
                 gradTiZ.data(), gradTiY.data(), cylsymm );
  std::cout << "thermal gradient interpolation gradTi " << gradTi[0] << " "
            << gradTi[1] << " " << gradTi[2] << " " << std::endl;

  // Applying background values at material boundaries
  std::for_each(boundaries.begin(), boundaries.end() - 1,
                boundary_init(background_Z, background_amu, nR_Dens, nZ_Dens,
                              DensGridr.data(), DensGridz.data(), ni.data(),
                              ne.data(), nR_Bfield, nZ_Bfield,
                              bfieldGridr.data(), bfieldGridz.data(), br.data(),
                              bz.data(), by.data(), nR_Temp, nZ_Temp,
                              TempGridr.data(), TempGridz.data(), ti.data(),
                              te.data(), biasPotential, biased_surface, surface_potential,
                              use_3d_geom, cylsymm ));

  std::cout << "Completed Boundary Init " << std::endl;
  std::cout << "periodicy "<<boundaries[nLines].periodic << std::endl;
  std::cout << "periodicx "<<boundaries[nLines].periodic_bc_x << std::endl;

  std::vector<double> surface_potential_out(nLines,0.0);
  std::vector<double> surface_te_out(nLines,0.0);
  std::vector<double> surface_ti_out(nLines,0.0);

  for( int i=0; i<nLines; i++)
{
  surface_potential_out[i] = boundaries[i].potential;
  surface_te_out[i] = boundaries[i].te;
  surface_ti_out[i] = boundaries[i].ti;
}
  netCDF::NcFile ncFile_boundary("output/boundary_values.nc", netCDF::NcFile::replace);
  netCDF::NcDim nc_nLine = ncFile_boundary.addDim("nBoundaries", nLines);

  netCDF::NcVar nc_bound_potential = ncFile_boundary.addVar("potential", netcdf_precision,nc_nLine);
  netCDF::NcVar nc_bound_te = ncFile_boundary.addVar("te", netcdf_precision,nc_nLine);
  netCDF::NcVar nc_bound_ti = ncFile_boundary.addVar("ti", netcdf_precision,nc_nLine);
  nc_bound_potential.putVar(&surface_potential_out[0]);
  nc_bound_te.putVar(&surface_te_out[0]);
  nc_bound_ti.putVar(&surface_ti_out[0]);
  ncFile_boundary.close();

  // Efield
  int nR_PreSheathEfield = 1;
  int nY_PreSheathEfield = 1;
  int nZ_PreSheathEfield = 1;
  int nPSEs = 1;
  std::string PSECfg = "backgroundPlasmaProfiles.Efield.";
// sim::Array<float> preSheathEGridy(1);

  sim::Array< gitr_precision > PSEr( nPSEs );
  sim::Array< gitr_precision > PSEz( nPSEs );
  sim::Array< gitr_precision > PSEt( nPSEs );

  sim::Array< gitr_precision > preSheathEGridr( nR_PreSheathEfield );
  sim::Array< gitr_precision > preSheathEGridy( nY_PreSheathEfield );
  sim::Array< gitr_precision > preSheathEGridz( nZ_PreSheathEfield );

  if( presheath_efield > 0 )
  {

  std::cout << "Using presheath Efield " << std::endl;
  if( presheath_interp == 1 )
  {
  nR_PreSheathEfield = nR_Lc;
  nY_PreSheathEfield = nY_Lc;
  nZ_PreSheathEfield = nZ_Lc;
  }

  std::string efieldFile;

  if( presheath_interp > 1 )
  {
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, PSECfg + "fileString", efieldFile);
    nR_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNrString");
    nZ_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNzString");

    if( presheath_interp > 2 )
    {
    nY_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNyString");
    }

#if USE_MPI > 0
  }
  MPI_Bcast(&nR_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  preSheathEGridr.resize( nR_PreSheathEfield );
  preSheathEGridy.resize( nY_PreSheathEfield );
  preSheathEGridz.resize( nZ_PreSheathEfield );

  nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

  PSEr.resize( nPSEs );
  PSEz.resize( nPSEs );
  PSEt.resize( nPSEs );

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    if( presheath_interp == 0 )
    {
    getVariable(cfg, PSECfg + "Er", PSEr[0]);
    getVariable(cfg, PSECfg + "Et", PSEt[0]);
    getVariable(cfg, PSECfg + "Ez", PSEz[0]);
    }
    else if( presheath_interp > 1 )
    {
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridRString",
                 preSheathEGridr[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridZString",
                 preSheathEGridz[0]);
  if( presheath_interp > 2 )
  {
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridYString",
                 preSheathEGridy[0]);
  }

  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "radialComponentString",
                 PSEr[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg,
                 "toroidalComponentString", PSEt[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "axialComponentString",
                 PSEz[0]);
    }

#if USE_MPI > 0
  }
  MPI_Bcast(preSheathEGridr.data(), nR_PreSheathEfield, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(preSheathEGridy.data(), nY_PreSheathEfield, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(preSheathEGridz.data(), nZ_PreSheathEfield, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(PSEr.data(), nPSEs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(PSEt.data(), nPSEs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(PSEz.data(), nPSEs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

if( presheath_interp == 1 )
{

  for (int i = 0; i < nR_PreSheathEfield; i++) {
    preSheathEGridr[i] = gridRLc[i];
    std::cout << "gridRLc " << gridRLc[i] << std::endl;
  }
  for (int i = 0; i < nY_PreSheathEfield; i++) {
    preSheathEGridy[i] = gridYLc[i];
  }
  for (int i = 0; i < nZ_PreSheathEfield; i++) {
    preSheathEGridz[i] = gridZLc[i];
  }
  std::cout << "length of PSE vec " << nPSEs << std::endl;
  gitr_precision teLocal1 = 0.0;
  gitr_precision BLocal1[3] = {0.0, 0.0, 0.0};
  gitr_precision Bnorm1[3] = {0.0, 0.0, 0.0};
  gitr_precision Bmag1 = 0.0;
  int index1 = 0;
  gitr_precision absS1 = 0.0;
  gitr_precision Epar = 0.0;
  for (int i = 0; i < nR_PreSheathEfield; i++) {
#if LC_INTERP == 3
    for (int k = 0; k < nY_PreSheathEfield; k++) {
      thisY = preSheathEGridy[k];
#endif
      for (int j = 0; j < nZ_PreSheathEfield; j++) {
        teLocal1 = interp2dCombined(preSheathEGridr[i], 0.0, preSheathEGridz[j],
                                    nR_Temp, nZ_Temp, &TempGridr.front(),
                                    &TempGridz.front(), &te.front(), cylsymm );
        interp2dVector(&BLocal1[0], gridRLc[i], 0.0, gridZLc[j], nR_Bfield,
                       nZ_Bfield, bfieldGridr.data(), bfieldGridz.data(),
                       br.data(), bz.data(), by.data(), cylsymm );
        Bmag1 = std::sqrt(BLocal1[0] * BLocal1[0] + BLocal1[1] * BLocal1[1] +
                     BLocal1[2] * BLocal1[2]);
        Bnorm1[0] = BLocal1[0] / Bmag1;
        Bnorm1[1] = BLocal1[1] / Bmag1;
        Bnorm1[2] = BLocal1[2] / Bmag1;

#if LC_INTERP == 3
        index1 = i + k * nR_PreSheathEfield +
                 j * nR_PreSheathEfield * nY_PreSheathEfield;
        // std::cout << "flowv calc index " << index << std::endl;
#else
      index1 = i + j * nR_PreSheathEfield;
#endif
        absS1 = std::abs(s[index1]);
        Epar = teLocal1 *
               (0.5 * Lc[index1] / absS1 /
                    std::sqrt(0.25 * Lc[index1] * Lc[index1] - absS1 * absS1) -
                1.0 / absS1);
        if (std::isnan(Epar))
          Epar = 0.0;
        PSEr[index1] = std::copysign(1.0,s[index1]) * Bnorm1[0] * Epar;
        PSEt[index1] = std::copysign(1.0,s[index1]) * Bnorm1[1] * Epar;
        PSEz[index1] = std::copysign(1.0,s[index1]) * Bnorm1[2] * Epar;
      }
#if LC_INTERP == 3
    }
#endif
  }
  sim::Array<gitr_precision> PSErSub(nPSEs), PSEzSub(nPSEs), PSEySub(nPSEs);

  for (int i = 0; i < nR_Lc; i++) {
    for (int j = 0; j < nY_Lc; j++) {
      for (int k = 0; k < nZ_Lc; k++) {
        int index = i + j * nR_Lc + k * nR_Lc * nY_Lc;
        if (noIntersectionNodes[index] == 1) {
          surroundingMinimumR = 0.0;
          surroundingMinimumY = 0.0;
          surroundingMinimumZ = 0.0;
          for (int ii = i - 1; ii < i + 2; ii++) {
            for (int jj = j - 1; jj < j + 2; jj++) {
              for (int kk = k - 1; kk < k + 2; kk++) {
                iterIndex = ii + jj * nR_Lc + kk * nR_Lc * nY_Lc;
                if (iterIndex > 0 && iterIndex < nFlowVs) {
                  if (noIntersectionNodes[iterIndex] == 0) {
                    if (std::abs(PSEr[iterIndex]) > std::abs(surroundingMinimumR)) {
                      surroundingMinimumR = PSEr[iterIndex];
                    }
                    if (std::abs(PSEt[iterIndex]) > std::abs(surroundingMinimumY)) {
                      surroundingMinimumY = PSEt[iterIndex];
                    }
                    if (std::abs(PSEz[iterIndex]) > std::abs(surroundingMinimumZ)) {
                      surroundingMinimumZ = PSEz[iterIndex];
                    }
                  }
                }
              }
            }
          }
          PSErSub[index] = surroundingMinimumR;
          PSEySub[index] = surroundingMinimumY;
          PSEzSub[index] = surroundingMinimumZ;
        }
      }
    }
  }
  for (int i = 0; i < nPSEs; i++) {
    if (i == 282839) {
      std::cout << " noIntersectionNodes " << noIntersectionNodes[i]
                << std::endl;
    }
    if (noIntersectionNodes[i] == 1) {
      PSEr[i] = PSErSub[i];
      PSEt[i] = PSEySub[i];
      PSEz[i] = PSEzSub[i];
    }
  }
  NcVar nc_PSEr = ncFileLC.addVar("PSEr", netcdf_precision, nc_nTracers);
  NcVar nc_PSEt = ncFileLC.addVar("PSEt", netcdf_precision, nc_nTracers);
  NcVar nc_PSEz = ncFileLC.addVar("PSEz", netcdf_precision, nc_nTracers);
  nc_PSEr.putVar(&PSEr[0]);
  nc_PSEt.putVar(&PSEt[0]);
  nc_PSEz.putVar(&PSEz[0]);
}
  }
  else
  {
  nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;
  sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
  sim::Array<gitr_precision> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);
  }

  std::string outnamePSEfieldR = "PSEfieldR.m";
  std::string outnamePSEfieldZ = "PSEfieldZ.m";
  std::string outnamePSEGridR = "PSEgridR.m";
  std::string outnamePSEGridZ = "PSEgridZ.m";

  std::cout << "Completed presheath Efield Init " << std::endl;
  sim::Array<gitr_precision> Efieldr(nR_Bfield * nZ_Bfield),
      Efieldz(nR_Bfield * nZ_Bfield), Efieldt(nR_Bfield * nZ_Bfield),
      minDist(nR_Bfield * nZ_Bfield);

  /* Captain! Likely dead block below */
  if( sheath_efield > 0 )
  {
    gitr_precision thisE0[3] = {0.0, 0.0, 0.0};
    gitr_precision minDist0 = 0.0;
    int minInd_bnd = 0;
    for (int i = 0; i < 1000; i++) {
        minDist0 =
            getE(0.0, 0.0, 1.0E-6*i, thisE0, boundaries.data(),
                nLines,
                nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
                n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
                &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
                &closeGeom_sheath.front(), minInd_bnd, biased_surface,
                use_3d_geom, geom_hash_sheath, cylsymm, sheath_model_type);
    }

    if( efield_interp == 2 )
      {
      int nR_dtsEfield, nZ_dtsEfield;

      int d1 = read_profileNs(
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNrString"),
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNzString"),
          nR_dtsEfield, nZ_dtsEfield);

      sim::Array<gitr_precision> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
      sim::Array<gitr_precision> dtsE(nR_dtsEfield * nZ_dtsEfield);

      int d2 = read_profile1d(
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridRString"),
          dtsEfieldGridr);

      int d3 = read_profile1d(
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridZString"),
          dtsEfieldGridz);

      int d4 = read_profile2d(
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
          cfg.lookup("backgroundPlasmaProfiles.dtsEfield.sheathDTS"), dtsE);
      }
  }
  else
  {
  int nR_dtsEfield = 1;
  int nZ_dtsEfield = 1;
  sim::Array<gitr_precision> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
  sim::Array<gitr_precision> dtsE(nR_dtsEfield * nZ_dtsEfield);
  }

  std::string outnameEfieldR = "EfieldR.m";
  std::string outnameEfieldZ = "EfieldZ.m";
  std::string outnameEfieldT = "EfieldT.m";
  std::string outnameMinDist = "DistToSurface.m";

  gitr_precision netX0 = 0.0;
  gitr_precision netX1 = 0.0;
  gitr_precision netY0 = 0.0;
  gitr_precision netY1 = 0.0;
  gitr_precision netZ0 = 0.0;
  gitr_precision netZ1 = 0.0;

  int net_nX = 0;
  int net_nY = 0;
  int net_nZ = 0;

  int nBins = 0;
  int nSpec = 0;

  size_t net_Bins_size = 0;
  size_t net_BinsTotal_size = 0;

  sim::Array< gitr_precision > gridX_bins(net_nX);
  sim::Array< gitr_precision > gridY_bins(net_nY);
  sim::Array< gitr_precision > gridZ_bins(net_nZ);

  if( spectroscopy > 0 )
  {
      if (world_rank == 0) {
        if (cfg.lookupValue("diagnostics.netx0", netX0) &&
            cfg.lookupValue("diagnostics.netx1", netX1) &&
            cfg.lookupValue("diagnostics.nety0", netY0) &&
            cfg.lookupValue("diagnostics.nety1", netY1) &&
            cfg.lookupValue("diagnostics.netz0", netZ0) &&
            cfg.lookupValue("diagnostics.netz1", netZ1) &&
            cfg.lookupValue("diagnostics.nX", net_nX) &&
            cfg.lookupValue("diagnostics.nY", net_nY) &&
            cfg.lookupValue("diagnostics.nZ", net_nZ) &&
            cfg.lookupValue("diagnostics.densityChargeBins", nBins)) {
          std::cout << "Spectroscopy net imported" << std::endl;
        } else {
          std::cout << "ERROR: Could not get spectroscopy net string info from "
                      "input file "
                    << std::endl;
        }
      }

      if( spectroscopy < 3 )
      {
        nSpec = (nBins + 1) * net_nX * net_nZ;
      }
      else
      {
        nSpec = (nBins + 1) * net_nX * net_nY * net_nZ;
      }

#if USE_MPI > 0
  MPI_Bcast(&netX0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netX1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netY0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netY1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netZ0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&netZ1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&net_nX, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&net_nY, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&net_nZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nBins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nSpec, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  std::cout << "spec bin Ns " << nBins << " " << net_nX << " " << net_nY << " "
            << net_nZ << std::endl;
  if( spectroscopy < 3 )
  {
  net_Bins_size = ((nBins + 1) * net_nX * net_nZ);
  net_BinsTotal_size = ((nBins + 1) * net_nX * net_nZ);
  }
  else
  {
  net_Bins_size = ((nBins + 1) * net_nX * net_nY * net_nZ);
  net_BinsTotal_size = ((nBins + 1) * net_nX * net_nY * net_nZ);
  }

  gridX_bins.resize(net_nX);
  gridY_bins.resize(net_nY);
  gridZ_bins.resize(net_nZ);

  for (int i = 0; i < net_nX; i++) {
    gridX_bins[i] = netX0 + 1.0 / (net_nX - 1) * i * (netX1 - netX0);
  }
  for (int i = 0; i < net_nY; i++) {
    gridY_bins[i] = netY0 + 1.0 / (net_nY - 1) * i * (netY1 - netY0);
  }

  for (int i = 0; i < net_nZ; i++) {
    gridZ_bins[i] = netZ0 + i * 1.0 / (net_nZ - 1) * (netZ1 - netZ0);
  }
}

  sim::Array<gitr_precision> net_Bins( net_Bins_size, 0.0 );
  sim::Array<gitr_precision> net_BinsTotal( net_BinsTotal_size, 0.0 );

  gitr_precision perpDiffusionCoeff = 0.0;
  if (world_rank == 0) {
    if (cfg.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                        perpDiffusionCoeff)) {
    } else {
      std::cout << "ERROR: could not get perpendicular diffusion coefficient "
                   "from input file"
                << std::endl;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&perpDiffusionCoeff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef __CUDACC__
  cout << "Using THRUST" << endl;
#else
  cout << "Not using THRUST" << endl;
#endif
  gitr_precision dt;
  int nP = 0;          // cfg.lookup("impurityParticleSource.nP");
  long nParticles = 0; // nP;
  int nT = 0;
  gitr_precision max_dt = 1.0e5;
  if (world_rank == 0) {
    nP = cfg.lookup("impurityParticleSource.nP");
    nParticles = nP;
    if (cfg.lookupValue("timeStep.dt", dt) &&
        cfg.lookupValue("timeStep.nT", nT)) {
      cout << "Number of time steps: " << nT << " With dt = " << dt << endl;
      cout << "Number of particles: " << nP << endl;
    } else {
      std::cout << "ERROR: could not get nT, dt, or nP from input file"
                << std::endl;
    }
    if (cfg.lookupValue("timeStep.max_dt", max_dt)) 
    {
      cout << "Max dt: " << max_dt << endl;
    } 
    else 
    {
      std::cout << "WARNING: maximum dt is not specified in input file, using default of 1.0e5 seconds"
                << std::endl;
    }
  }

#if USE_MPI > 0
  MPI_Bcast(&dt, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nP, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nParticles, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  sim::Array<int> nPPerRank(world_size, 0), pStartIndx(world_size, 0),
      pDisplacement(world_size, 0), pHistPerNode(world_size, 0),
      nActiveParticlesOnRank(world_size, 0);
  int countP = 0;
  if (nP >= world_size) {
    for (int i = 0; i < world_size; i++) {
      nPPerRank[i] = std::floor(nP / world_size);
      if (i == 0) {
        nPPerRank[i] = nPPerRank[i] + nP % world_size;
      }
      pStartIndx[i] = countP;
      countP = countP + nPPerRank[i];
      std::cout << "countP " << countP << std::endl;
    }
  } else {
    for (int i = 0; i < nP; i++) {
      nPPerRank[i] = 1;
      pStartIndx[i] = countP;
      countP = countP + nPPerRank[i];
    }
  }

  for (int i = 0; i < world_size; i++) {
    nActiveParticlesOnRank[i] = nPPerRank[i];
  }
         
  auto particleArray = new Particles(nParticles,1,cfg,gitr_flags);
  std::cout << "creating particle source" << std::endl;
  nP = cfg.lookup("impurityParticleSource.nP");
  sim::Array<gitr_precision>  px(nP), py(nP), pz(nP), pvx(nP), pvy(nP), pvz(nP), pamu(nP), pZ(nP), pcharge(nP);
  sim::Array<int> speciesType(nP);
  
  libconfig::Setting& speciesArray = cfg.lookup("impurityParticleSource.initialConditions.species");
  int num_species = speciesArray.getLength();
  printf("Number of species: %d\n",num_species);

  //  read in particle source or generate source 
  ParticleData particleData;

  if (particle_source_file == 1) {
      std::string ncParticleSourceFile;
      getVariable(cfg, "particleSource.ncFileString", ncParticleSourceFile);
      particleData = readParticleData(ncParticleSourceFile);
  } else if (particle_source_file == 0) {
      particleData = generateParticleData(cfg, particle_source_energy);
  } else {
      std::cerr << "Error: Invalid value for particle_source_file." << std::endl;
  }
  initializeParticleArray(particleData, particleArray, px, py, pz, pvx, pvy, pvz, pZ, pamu, pcharge, dt, speciesType);

  std::cout << "writing particles out file" << std::endl;
  netCDF::NcFile ncFile_particles("output/particleSource.nc", netCDF::NcFile::replace);
  netCDF::NcDim pNP = ncFile_particles.addDim("nP", nP);
  netCDF::NcVar p_surfNormx = ncFile_particles.addVar("surfNormX", netcdf_precision, pNP);
  netCDF::NcVar p_surfNormy = ncFile_particles.addVar("surfNormY", netcdf_precision, pNP);
  netCDF::NcVar p_surfNormz = ncFile_particles.addVar("surfNormZ", netcdf_precision, pNP);
  netCDF::NcVar p_vx = ncFile_particles.addVar("vx", netcdf_precision, pNP);
  netCDF::NcVar p_vy = ncFile_particles.addVar("vy", netcdf_precision, pNP);
  netCDF::NcVar p_vz = ncFile_particles.addVar("vz", netcdf_precision, pNP);
  netCDF::NcVar p_x = ncFile_particles.addVar("x", netcdf_precision, pNP);
  netCDF::NcVar p_y = ncFile_particles.addVar("y", netcdf_precision, pNP);
  netCDF::NcVar p_z = ncFile_particles.addVar("z", netcdf_precision, pNP);
  netCDF::NcVar p_charge = ncFile_particles.addVar("charge", netcdf_precision, pNP);
  netCDF::NcVar p_amu = ncFile_particles.addVar("amu", netcdf_precision, pNP);
  netCDF::NcVar p_Z = ncFile_particles.addVar("Z", netcdf_precision, pNP);
  // add species type
  netCDF::NcVar p_speciesType = ncFile_particles.addVar("speciesType", netCDF::ncInt, pNP);
  p_vx.putVar(&pvx[0]);
  p_vy.putVar(&pvy[0]);
  p_vz.putVar(&pvz[0]);
  p_x.putVar(&px[0]);
  p_y.putVar(&py[0]);
  p_z.putVar(&pz[0]);
  p_charge.putVar(&pcharge[0]);
  p_amu.putVar(&pamu[0]);
  p_Z.putVar(&pZ[0]);
  p_speciesType.putVar(&speciesType[0]);
  ncFile_particles.close();

  std::cout << "finished writing particles out file" << std::endl;
  int subSampleFac = 1;
  if (world_rank == 0) {
    if (cfg.lookupValue("diagnostics.trackSubSampleFactor", subSampleFac)) {
      std::cout << "Tracks subsample factor imported " << subSampleFac
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get tracks sub sample info from input file "
          << std::endl;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&subSampleFac, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  int nHistoriesPerParticle = (nT / subSampleFac) + 1;
  int nHistories = nHistoriesPerParticle * nP;
  sim::Array<gitr_precision> positionHistoryX(nHistories);
  sim::Array<gitr_precision> positionHistoryXgather(nHistories, 0.0);
  sim::Array<gitr_precision> positionHistoryY(nHistories);
  sim::Array<gitr_precision> positionHistoryYgather(nHistories);
  sim::Array<gitr_precision> positionHistoryZ(nHistories);
  sim::Array<gitr_precision> positionHistoryZgather(nHistories);
  sim::Array<gitr_precision> velocityHistory(nHistories);
  sim::Array<gitr_precision> velocityHistoryX(nHistories);
  sim::Array<gitr_precision> velocityHistoryY(nHistories);
  sim::Array<gitr_precision> velocityHistoryZ(nHistories);
  sim::Array<gitr_precision> velocityHistorygather(nHistories);
  sim::Array<gitr_precision> velocityHistoryXgather(nHistories);
  sim::Array<gitr_precision> velocityHistoryYgather(nHistories);
  sim::Array<gitr_precision> velocityHistoryZgather(nHistories);
  sim::Array<gitr_precision> chargeHistory(nHistories);
  sim::Array<gitr_precision> chargeHistoryGather(nHistories);
    sim::Array<gitr_precision> ZHistory(nHistories);
  sim::Array<gitr_precision> ZHistoryGather(nHistories);
  sim::Array<gitr_precision> weightHistory(nHistories);
  sim::Array<gitr_precision> weightHistoryGather(nHistories);

  if( particle_tracks > 0 )
  {
    for (int i = 0; i < world_size; i++) {
      pDisplacement[i] = pStartIndx[i] * nHistoriesPerParticle;
      pHistPerNode[i] = nPPerRank[i] * nHistoriesPerParticle;
    }
  const int *displ = &pDisplacement[0];
  const int *phpn = &pHistPerNode[0];
  std::cout << "History array length " << nHistories << std::endl;
  }

  gitr_precision *finalPosX = new gitr_precision[nP];
  gitr_precision *finalPosY = new gitr_precision[nP];
  gitr_precision *finalPosZ = new gitr_precision[nP];
  gitr_precision *finalVx = new gitr_precision[nP];
  gitr_precision *finalVy = new gitr_precision[nP];
  gitr_precision *finalVz = new gitr_precision[nP];
  gitr_precision *transitTime = new gitr_precision[nP];
  gitr_precision *hitWall = new gitr_precision[nP];

  std::cout << "Beginning random number seeds" << std::endl;
  std::uniform_real_distribution<gitr_precision> dist(0, 1e6);
  thrust::counting_iterator<std::size_t> particleBegin(pStartIndx[world_rank]);
  thrust::counting_iterator<std::size_t> particleEnd(pStartIndx[world_rank] + nActiveParticlesOnRank[world_rank] );
  thrust::counting_iterator<std::size_t> particleOne(1);
  thrust::counting_iterator<std::size_t> particleZero(0);
  auto randInitStart_clock = gitr_time::now();

#ifdef __CUDACC__
     typedef curandState rand_type;
#else
     typedef std::mt19937 rand_type;
#endif
#if USE_CUDA
  sim::Array<rand_type> state1(nParticles);
#else
  sim::Array<rand_type> state1(nParticles);
#endif
    if( ionization > 0 ||
        perp_diffusion > 0 ||
        coulomb_collisions > 0 ||
        surface_model > 0 )
    {
#if USE_CUDA

  int *dev_i;
  cudaMallocManaged(&dev_i, sizeof(int));
  dev_i[0] = 0;
  std::cout << " about to do curandInit" << std::endl;
  thrust::for_each(thrust::device, particleBegin, particleEnd,
                   curandInitialize<>( &state1.front(), fixed_seeds ));
  std::cout << " finished curandInit" << std::endl;
#else
  std::random_device randDeviceInit;

  thrust::for_each( thrust::device, particleBegin, particleEnd,
                    curandInitialize<>( &state1.front(), fixed_seeds ) );

#endif
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
}
  auto randInitEnd_clock = gitr_time::now();
  std::chrono::duration<gitr_precision> fsRandInit = randInitEnd_clock - randInitStart_clock;
  printf(
      "Random Number Initialize time for node %i          is %6.3f (secs) \n",
      world_rank, fsRandInit.count());

  gitr_precision moveTime = 0.0;
  gitr_precision geomCheckTime = 0.0;
  gitr_precision ionizTime = 0.0;
#if USE_CUDA > 0
  int *dev_tt;
  cudaMallocManaged(&dev_tt, sizeof(int));
#else
  int *dev_tt = new int[1];
  *dev_tt = 0;
#endif
  int tt = 0;

  move_boris move_boris0(
      particleArray, dt, boundaries.data(), nLines, nR_Bfield, nZ_Bfield,
      bfieldGridr.data(), &bfieldGridz.front(), &br.front(), &bz.front(),
      &by.front(), nR_PreSheathEfield, nY_PreSheathEfield, nZ_PreSheathEfield,
      &preSheathEGridr.front(), &preSheathEGridy.front(),
      &preSheathEGridz.front(), &PSEr.front(), &PSEz.front(), &PSEt.front(),
      nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
      n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
      &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
      &closeGeom_sheath.front(),gitr_flags, sheath_efield, presheath_efield, biased_surface,
      geom_hash_sheath,
      use_3d_geom,
      cylsymm,
      max_dt, sheath_model_type);

  geometry_check geometry_check0(
      particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
      nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
      n_closeGeomElements.data(), &closeGeomGridr.front(),
      &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
      nEdist, E0dist, Edist, nAdist, A0dist, Adist, flux_ea, surface_model,
      geom_hash,
      use_3d_geom,
      cylsymm );

  sortParticles sort0(particleArray, nP,dev_tt, 10000,
                      nActiveParticlesOnRank.data());
  spec_bin spec_bin0(gitr_flags,particleArray, nBins, net_nX, net_nY, net_nZ,
                     &gridX_bins.front(), &gridY_bins.front(),
                     &gridZ_bins.front(), &net_Bins.front(), dt, cylsymm, spectroscopy );

  gitr_precision *uni;

  if( ionization > 0 )
  {
#if USE_CUDA > 0
  cudaMallocManaged(&uni, sizeof(gitr_precision));
#else
  uni = new gitr_precision[1];
  *uni = 0;
#endif
  }

  ionize<rand_type> ionize0(
      gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(),
      uni, cylsymm );

  recombine<rand_type> recombine0(
      particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
       &TempGridr.front(),
      &TempGridz.front(), &te.front(),
      gitr_flags, cylsymm );

  crossFieldDiffusion crossFieldDiffusion0( gitr_flags,
      particleArray, dt, &state1.front(), perpDiffusionCoeff, nR_Bfield,
      nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(),
      &bz.front(), &by.front(), perp_diffusion, cylsymm );

  coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, bfieldGridr.data(),
      &bfieldGridz.front(), &br.front(), &bz.front(), &by.front(),gitr_flags, flowv_interp,
      cylsymm, field_aligned_values );

  thermalForce thermalForce0(gitr_flags,
      particleArray, dt, background_amu, nR_gradT, nZ_gradT, gradTGridr.data(),
      gradTGridz.data(), gradTiR.data(), gradTiZ.data(), gradTiY.data(),
      gradTeR.data(), gradTeZ.data(), gradTeY.data(), nR_Bfield, nZ_Bfield,
      bfieldGridr.data(), &bfieldGridz.front(), &br.front(), &bz.front(),
      &by.front(), cylsymm );

  reflection reflection0( particleArray, dt, &state1.front(), nLines, &boundaries[0], surfaces, flux_ea, use_3d_geom, cylsymm, nspecies);

  history history0(particleArray, dev_tt, nT, subSampleFac, nP,
                   &positionHistoryX.front(), &positionHistoryY.front(),
                   &positionHistoryZ.front(), &velocityHistory.front(),
                   &velocityHistoryX.front(), &velocityHistoryY.front(),
                   &velocityHistoryZ.front(), &chargeHistory.front(), &ZHistory.front(),
                   &weightHistory.front());


// if( force_eval > 0 )
//   {
//   if (world_rank == 0) {
//     int nR_force, nZ_force;
//     gitr_precision forceX0, forceX1, forceZ0, forceZ1, testEnergy;
//     std::string forceCfg = "forceEvaluation.";

//     getVariable(cfg, forceCfg + "nR", nR_force);
//     getVariable(cfg, forceCfg + "nZ", nZ_force);
//     std::vector<gitr_precision> forceR(nR_force, 0.0), forceZ(nZ_force, 0.0);
//     std::vector<gitr_precision> tIon(nR_force * nZ_force, 0.0),
//         tRecomb(nR_force * nZ_force, 0.0);
//     std::vector<gitr_precision> dvEr(nR_force * nZ_force, 0.0),
//         dvEz(nR_force * nZ_force, 0.0), dvEt(nR_force * nZ_force, 0.0);
//     std::vector<gitr_precision> dvBr(nR_force * nZ_force, 0.0),
//         dvBz(nR_force * nZ_force, 0.0), dvBt(nR_force * nZ_force, 0.0);
//     std::vector<gitr_precision> dvCollr(nR_force * nZ_force, 0.0),
//         dvCollz(nR_force * nZ_force, 0.0), dvCollt(nR_force * nZ_force, 0.0);
//     std::vector<gitr_precision> dvITGr(nR_force * nZ_force, 0.0),
//         dvITGz(nR_force * nZ_force, 0.0), dvITGt(nR_force * nZ_force, 0.0);
//     std::vector<gitr_precision> dvETGr(nR_force * nZ_force, 0.0),
//         dvETGz(nR_force * nZ_force, 0.0), dvETGt(nR_force * nZ_force, 0.0);
//     getVariable(cfg, forceCfg + "X0", forceX0);
//     getVariable(cfg, forceCfg + "X1", forceX1);
//     getVariable(cfg, forceCfg + "Z0", forceZ0);
//     getVariable(cfg, forceCfg + "Z1", forceZ1);
//     getVariable(cfg, forceCfg + "particleEnergy", testEnergy);
//     for (int i = 0; i < nR_force; i++) {
//       forceR[i] = forceX0 + (forceX1 - forceX0) * i / (nR_force - 1);
//     }
//     for (int i = 0; i < nZ_force; i++) {
//       forceZ[i] = forceZ0 + (forceZ1 - forceZ0) * i / (nZ_force - 1);
//     }
//     gitr_precision Btotal = 0.0;
//     for (int i = 0; i < nR_force; i++) {
//       for (int j = 0; j < nZ_force; j++) {
//         interp2dVector(&Btest[0], forceR[i], 0.0, forceZ[j], nR_Bfield,
//                        nZ_Bfield, bfieldGridr.data(), bfieldGridz.data(),
//                        br.data(), bz.data(), by.data(), cylsymm );
//         Btotal = vectorNorm(Btest);
//         gitr_precision testTi =
//             interp2dCombined(0.0, 0.1, 0.0, nR_Temp, nZ_Temp, TempGridr.data(),
//                              TempGridz.data(), ti.data(), cylsymm );
//         particleArray->setParticle(0, forceR[i], 0.0, forceZ[j], testTi, 0.0,
//                                    0.0, pZ[0], pamu[0], pcharge[0] + 1.0,1);
//         move_boris0(0);

//         if( ionization > 0 )
//         {
//           thrust::for_each(thrust::device,particleBegin,particleBegin,ionize0);
// 	        thrust::for_each(thrust::device,particleBegin,particleBegin,recombine0);
//         }

//         if( coulomb_collisions > 0 )
//         {
//         thrust::for_each(thrust::device,particleBegin,particleBegin,coulombCollisions0);
//         }

//         if( thermal_force > 0 )
//         {
//         thrust::for_each(thrust::device,particleBegin,particleBegin,thermalForce0);
//         }
//         dvEr[j * nR_force + i] = move_boris0.electricForce[0];
//         dvEz[j * nR_force + i] = move_boris0.electricForce[2];
//         dvEt[j * nR_force + i] = move_boris0.electricForce[1];
//         dvBr[j * nR_force + i] = move_boris0.magneticForce[0];
//         dvBz[j * nR_force + i] = move_boris0.magneticForce[2];
//         dvBt[j * nR_force + i] = move_boris0.magneticForce[1];

//         if( ionization > 0 )
//         {
//           tIon[j * nR_force + i] = ionize0.tion;
//           tRecomb[j * nR_force + i] = recombine0.tion;
//         }

//         if( coulomb_collisions > 0 )
//         {
//         dvCollr[j * nR_force + i] = coulombCollisions0.dv[0];
//         dvCollz[j * nR_force + i] = coulombCollisions0.dv[2];
//         dvCollt[j * nR_force + i] = coulombCollisions0.dv[1];
//         }
//         if( thermal_force > 0 )
//         {
//         dvITGr[j * nR_force + i] = thermalForce0.dv_ITGx;
//         dvITGz[j * nR_force + i] = thermalForce0.dv_ITGz;
//         dvITGt[j * nR_force + i] = thermalForce0.dv_ITGy;
//         dvETGr[j * nR_force + i] = thermalForce0.dv_ETGx;
//         dvETGz[j * nR_force + i] = thermalForce0.dv_ETGz;
//         dvETGt[j * nR_force + i] = thermalForce0.dv_ETGy;
//         }
//       }
//     }
//     std::cout << " about to write ncFile_forces " << std::endl;
//     netCDF::NcFile ncFile_force("output/forces.nc", netCDF::NcFile::replace);
//     netCDF::NcDim nc_nRf = ncFile_force.addDim("nR", nR_force);
//     netCDF::NcDim nc_nZf = ncFile_force.addDim("nZ", nZ_force);
//     vector<netCDF::NcDim> forceDims;
//     forceDims.push_back(nc_nZf);
//     forceDims.push_back(nc_nRf);
//     netCDF::NcVar forceRf = ncFile_force.addVar("r", netcdf_precision, nc_nRf);
//     netCDF::NcVar forceZf = ncFile_force.addVar("z", netcdf_precision, nc_nZf);
//     netCDF::NcVar nction = ncFile_force.addVar("tIon", netcdf_precision, forceDims);
//     netCDF::NcVar nctrec = ncFile_force.addVar("tRec", netcdf_precision, forceDims);
//     netCDF::NcVar dvErf = ncFile_force.addVar("dvEr", netcdf_precision, forceDims);
//     netCDF::NcVar dvEzf = ncFile_force.addVar("dvEz", netcdf_precision, forceDims);
//     netCDF::NcVar dvEtf = ncFile_force.addVar("dvEt", netcdf_precision, forceDims);
//     netCDF::NcVar dvBrf = ncFile_force.addVar("dvBr", netcdf_precision, forceDims);
//     netCDF::NcVar dvBzf = ncFile_force.addVar("dvBz", netcdf_precision, forceDims);
//     netCDF::NcVar dvBtf = ncFile_force.addVar("dvBt", netcdf_precision, forceDims);
//     netCDF::NcVar dvCollrf = ncFile_force.addVar("dvCollr", netcdf_precision, forceDims);
//     netCDF::NcVar dvCollzf = ncFile_force.addVar("dvCollz", netcdf_precision, forceDims);
//     netCDF::NcVar dvColltf = ncFile_force.addVar("dvCollt", netcdf_precision, forceDims);
//     netCDF::NcVar dvITGrf = ncFile_force.addVar("dvITGr", netcdf_precision, forceDims);
//     netCDF::NcVar dvITGzf = ncFile_force.addVar("dvITGz", netcdf_precision, forceDims);
//     netCDF::NcVar dvITGtf = ncFile_force.addVar("dvITGt", netcdf_precision, forceDims);
//     netCDF::NcVar dvETGrf = ncFile_force.addVar("dvETGr", netcdf_precision, forceDims);
//     netCDF::NcVar dvETGzf = ncFile_force.addVar("dvETGz", netcdf_precision, forceDims);
//     netCDF::NcVar dvETGtf = ncFile_force.addVar("dvETGt", netcdf_precision, forceDims);
//     forceRf.putVar(&forceR[0]);
//     forceZf.putVar(&forceZ[0]);
//     nction.putVar(&tIon[0]);
//     nctrec.putVar(&tRecomb[0]);
//     dvErf.putVar(&dvEr[0]);
//     dvEzf.putVar(&dvEz[0]);
//     dvEtf.putVar(&dvEt[0]);
//     dvBrf.putVar(&dvBr[0]);
//     dvBzf.putVar(&dvBz[0]);
//     dvBtf.putVar(&dvBt[0]);
//     dvCollrf.putVar(&dvCollr[0]);
//     dvCollzf.putVar(&dvCollz[0]);
//     dvColltf.putVar(&dvCollt[0]);
//     dvITGrf.putVar(&dvITGr[0]);
//     dvITGzf.putVar(&dvITGz[0]);
//     dvITGtf.putVar(&dvITGt[0]);
//     dvETGrf.putVar(&dvETGr[0]);
//     dvETGzf.putVar(&dvETGz[0]);
//     dvETGtf.putVar(&dvETGt[0]);
//     ncFile_force.close();
//     particleArray->setParticleV(0, px[0], py[0], pz[0], pvx[0], pvy[0], pvz[0],
//                                 pZ[0], pamu[0], pcharge[0], dt);
//   }
// }


  auto start_clock = gitr_time::now();
  std::chrono::duration<gitr_precision> fs1 = start_clock - gitr_start_clock;
  printf("Initialize time for node %i          is %6.3f (secs) \n", world_rank,
         fs1.count());
  gitr_precision testFlowVec[3] = {0.0};
  if( field_aligned_values > 0 )
  {
  interpFieldAlignedVector(&testFlowVec[0], 1.4981, 0.0, 1.0, nR_flowV,
                           nZ_flowV, flowVGridr.data(), flowVGridz.data(),
                           flowVr.data(), flowVz.data(), flowVt.data(),
                           nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                           bfieldGridz.data(), br.data(), bz.data(), by.data(), cylsymm );
  }
  else
  {
  interp2dVector(&testFlowVec[0], 1.4981, 0.0, 1.0, nR_flowV, nZ_flowV,
                 flowVGridr.data(), flowVGridz.data(), flowVr.data(),
                 flowVz.data(), flowVt.data(), cylsymm );
  }

  gitr_precision leakZ = 0.0;
  if (world_rank == 0) {

    std::string diagnosticCfg = "diagnostics.";

    getVariable(cfg, diagnosticCfg + "leakZ", leakZ);
  }
#if USE_MPI > 0
  MPI_Bcast(&leakZ, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  for (int i = 0; i < nP; i++)
    particleArray->leakZ[i] = leakZ;

  std::cout << "Flow velocities " << testFlowVec[0] << " " << testFlowVec[1]
            << " " << testFlowVec[2] << std::endl;
 

  std::cout << "Starting main loop" << std::endl;
  // Main time loop
#if __CUDACC__
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    sim::Array<int> tmpInt(1, 1), tmpInt2(1, 1);
#ifdef __CUDACC__
    cudaDeviceSynchronize();
#endif
    for (tt; tt < nT; tt++) {

       
       if( sort_particles > 0 )
       {
       dev_tt[0] = tt;
      //std::cout << " tt " << tt << std::endl;
      thrust::for_each(thrust::host, particleBegin, particleOne, sort0);
      particleEnd = particleZero + nActiveParticlesOnRank[0];
#ifdef __CUDACC__
      cudaDeviceSynchronize();
#endif
       }

      if( particle_tracks > 0 )
      {
      thrust::for_each(thrust::device, particleBegin, particleEnd, history0);
      }

      thrust::for_each(thrust::device, particleBegin, particleEnd, move_boris0);

#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       geometry_check0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif

      if( spectroscopy > 0 )
      {
      thrust::for_each(thrust::device, particleBegin, particleEnd, spec_bin0);
      }
      if( ionization > 0 )
      {
        thrust::for_each(thrust::device, particleBegin, particleEnd,ionize0);
        thrust::for_each(thrust::device, particleBegin, particleEnd, recombine0);
      }

      if( perp_diffusion > 0 )
      {
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       crossFieldDiffusion0);
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       geometry_check0);
      }

      if( coulomb_collisions > 0 )
      {
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       coulombCollisions0);
      }

      if( thermal_force > 0 )
      {
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       thermalForce0);
      }

      if( surface_model > 0 )
      {
          thrust::for_each(thrust::device, particleBegin, particleEnd, reflection0);
      }
    }
  
    if( particle_tracks > 0 )
    {
    tt = nT;
    thrust::for_each(thrust::device, particleBegin, particleEnd, history0);
    }


  // Ensure that all time step loop GPU kernels are complete before proceeding
#ifdef __CUDACC__
  cudaDeviceSynchronize();
#endif
  auto finish_clock = gitr_time::now();
  std::chrono::duration<gitr_precision> fs = finish_clock - start_clock;
  printf("Time taken          is %6.3f (secs) \n", fs.count());
  printf("Time taken per step is %6.3f (secs) \n", fs.count() / (gitr_precision)nT);

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

  // Get the number of unique species - each nuclear charge is a species
  std::map<double, std::string> uniqueSpecies;
  for (int i = 0; i < nP; i++) {
      double charges = particleArray->Z[i];
      std::string speciesNames = materialData[particleArray->Z[i]].name; 
      uniqueSpecies[charges] = speciesNames; 
  }

  // Display number of unique species
  std::cout << "Number of unique species: " << uniqueSpecies.size() << std::endl;

  // // Iterate over the unique species and save the data
  // for (const auto& pair : uniqueSpecies) {
  //   std::cout << "Unique charge: " << pair.first << ", Species Name: " << pair.second << std::endl;

  //   // Assuming you want the filename to be based on the species name
  //   std::string filename = "output/surface." + pair.second + ".nc";

  //   // Create and open a new NetCDF file
  //   netCDF::NcFile ncFile(filename, netCDF::NcFile::replace);

  //   int nSpecies = uniqueSpecies.size();

  //   // Create dimensions
  //   netCDF::NcDim nc_nLines = ncFile.addDim("nSurfaces", nSurfaces);
  //   netCDF::NcDim nc_nSpecies = ncFile.addDim("nSpecies", nSpecies);
  //   netCDF::NcDim nc_nEnergies = ncFile.addDim("nEnergies", nEdist);
  //   netCDF::NcDim nc_nAngles = ncFile.addDim("nAngles", nAdist);

  //   // Create and define variables
  //   std::vector<netCDF::NcDim> dims1 = {nc_nLines};
  //   std::vector<netCDF::NcDim> dimsSurfE = {nc_nLines, nc_nSpecies, nc_nEnergies, nc_nAngles};

  //   // Define variables
  //   netCDF::NcVar nc_grossDep = ncFile.addVar("grossDeposition", netcdf_precision, dims1);
  //   netCDF::NcVar nc_grossEro = ncFile.addVar("grossErosion", netcdf_precision, dimsSurfE);
  //   netCDF::NcVar nc_aveSpyl = ncFile.addVar("aveSpyl", netcdf_precision, dims1);
  //   netCDF::NcVar nc_spylCounts = ncFile.addVar("spylCounts", netCDF::ncInt, dims1);
  //   netCDF::NcVar nc_surfNum = ncFile.addVar("surfaceNumber", netCDF::ncInt, dims1);
  //   netCDF::NcVar nc_sumParticlesStrike = ncFile.addVar("sumParticlesStrike", netCDF::ncInt, dims1);
  //   netCDF::NcVar nc_sumWeightStrike = ncFile.addVar("sumWeightStrike", netcdf_precision, dims1);
  //   netCDF::NcVar nc_surfEDist = ncFile.addVar("surfEDist", netcdf_precision, dimsSurfE);
  //   netCDF::NcVar nc_surfReflDist = ncFile.addVar("surfReflDist", netcdf_precision, dimsSurfE);
  //   netCDF::NcVar nc_surfSputtDist = ncFile.addVar("surfSputtDist", netcdf_precision, dimsSurfE);

  //   // Collect and manage the data
  //   std::vector<int> surfaceNumbers(nSurfaces, 0);
  //   std::vector<std::vector<gitr_precision>> multiSpeciesGrossErosion(nSurfaces, std::vector<gitr_precision>(nSpecies, 0));

  //   int srf = 0;
  //   for (int i = 0; i < nLines; i++) {
  //       if (boundaries[i].surface) {
  //           surfaceNumbers[srf] = i;
  //           for (int j = 0; j < nSpecies; j++) {
  //               multiSpeciesGrossErosion[srf][j] = surfaces->grossErosion[srf] + multiSpeciesGrossErosion[srf][j];


  //           }
  //           srf++;
  //       }
  //   }  

  //   // Write data to NetCDF file
  //   nc_grossDep.putVar(&surfaces->grossDeposition[0]);
  //   nc_surfNum.putVar(&surfaceNumbers[0]);
  //   nc_grossEro.putVar(&multiSpeciesGrossErosion[0][0]);
  //   nc_aveSpyl.putVar(&surfaces->aveSputtYld[0]);
  //   nc_spylCounts.putVar(&surfaces->sputtYldCount[0]);
  //   nc_sumParticlesStrike.putVar(&surfaces->sumParticlesStrike[0]);
  //   nc_sumWeightStrike.putVar(&surfaces->sumWeightStrike[0]);
  //   nc_surfEDist.putVar(&surfaces->energyDistribution[0]);
  //   nc_surfReflDist.putVar(&surfaces->reflDistribution[0]);
  //   nc_surfSputtDist.putVar(&surfaces->sputtDistribution[0]);

  //   // Close the file
  //   ncFile.close();
  // }
////////////////////////////////////
////////////////////////////////////
////////////////////////////////////
  // float tmp202 =0.0;
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
// show memory usage of GPU
#if USE_CUDA
  if (world_rank == 0) {
    size_t free_byte;
    size_t total_byte;
    cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

    if (cudaSuccess != cuda_status) {

      printf("Error: cudaMemGetInfo fails, %s \n",
             cudaGetErrorString(cuda_status));
      exit(1);
    }

    double free_db = (double)free_byte;
    double total_db = (double)total_byte;
    double used_db = total_db - free_db;

    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
           used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0,
           total_db / 1024.0 / 1024.0);
    int nDevices;
    int nThreads;
    cudaGetDeviceCount(&nDevices);
    std::cout << "number of devices gotten " << nDevices << std::endl;
    for (int i = 0; i < nDevices; i++) {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("Device Number: %d\n", i);
      printf("  Device name: %s\n", prop.name);
      printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
      printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
      printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
             2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
      printf("  Total number of threads: %d\n",
             prop.maxThreadsPerMultiProcessor);
      nThreads = prop.maxThreadsPerMultiProcessor;
    }
  }
#endif
  std::cout << "reached gather initialization " << nP << std::endl;
  sim::Array<gitr_precision> xGather(nP, 0.0);
  sim::Array<gitr_precision> test0Gather(nP, 0.0);
  sim::Array<gitr_precision> test1Gather(nP, 0.0);
  sim::Array<gitr_precision> yGather(nP, 0.0);
  sim::Array<gitr_precision> zGather(nP, 0.0);
  sim::Array<gitr_precision> vGather(nP, 0.0);
  sim::Array<gitr_precision> vxGather(nP, 0.0);
  sim::Array<gitr_precision> vyGather(nP, 0.0);
  sim::Array<gitr_precision> vzGather(nP, 0.0);
  sim::Array<gitr_precision> hitWallGather(nP, 0.0);
  sim::Array<int> surfaceHitGather(nP, 0.0);
  sim::Array<gitr_precision> weightGather(nP, 0.0);
  sim::Array<gitr_precision> chargeGather(nP, 0.0);
  sim::Array<gitr_precision> firstIonizationTGather(nP, 0.0);
  sim::Array<gitr_precision> firstIonizationZGather(nP, 0.0);
  sim::Array<int> hasLeakedGather(nP, 0);

  std::cout << "Reached MPI barrier for gather" << std::endl;
  std::cout << "gather pstart and npperrank " << pStartIndx[world_rank] << " " << nPPerRank[world_rank] << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&particleArray->x[pStartIndx[world_rank]],nPPerRank[world_rank],
             MPI_FLOAT, &xGather[0], nPPerRank[world_rank], MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  std::cout << "Passed x"<< world_rank << std::endl;
  MPI_Gather(&particleArray->y[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &yGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->z[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &zGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->v[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &vGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->vx[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &vxGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->vy[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &vyGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->vz[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &vzGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&particleArray->hitWall[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &hitWallGather[0], nP / world_size,
             MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->hitWall[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &surfaceHitGather[0], nP / world_size,
             MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->weight[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &weightGather[0], nP / world_size,
             MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->charge[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &chargeGather[0], nP / world_size,
             MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->hasLeaked[world_rank * nP / world_size],
             nP / world_size, MPI_INT, &hasLeakedGather[0], nP / world_size,
             MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->firstIonizationT[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &firstIonizationTGather[0],
             nP / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->firstIonizationZ[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &firstIonizationZGather[0],
             nP / world_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->test0[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &test0Gather[0], nP / world_size,
             MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gather(&particleArray->test1[world_rank * nP / world_size],
             nP / world_size, MPI_FLOAT, &test1Gather[0], nP / world_size,
             MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if( particle_tracks > 0 )
  {

  std::vector<gitr_precision> exampleArray(4, 0.0);
  std::vector<gitr_precision> exampleArrayGather(4, 0.0);
  if (world_rank == 0) {
    exampleArray[0] = 1;
    exampleArray[1] = 1;
  }
  if (world_rank == 1) {
    exampleArray[2] = 2;
    exampleArray[3] = 2;
  }
  std::vector<int> exCount(2, 2), exDispl(2, 0);
  exDispl[0] = 0;
  exDispl[1] = 2;
  const int *exdispl = &exDispl[0];
  const int *excount = &exCount[0];

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gatherv(&positionHistoryX[pDisplacement[world_rank]],
              pHistPerNode[world_rank], MPI_FLOAT, &positionHistoryXgather[0],
              phpn, displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&positionHistoryY[pDisplacement[world_rank]],
              pHistPerNode[world_rank], MPI_FLOAT, &positionHistoryYgather[0],
              phpn, displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&positionHistoryZ[pDisplacement[world_rank]],
              pHistPerNode[world_rank], MPI_FLOAT, &positionHistoryZgather[0],
              phpn, displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&velocityHistory[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &velocityHistorygather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(&velocityHistoryX[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &velocityHistoryXgather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(&velocityHistoryY[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &velocityHistoryYgather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(&velocityHistoryZ[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &velocityHistoryZgather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(&chargeHistory[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &chargeHistoryGather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Gatherv(&weightHistory[pStartIndx[world_rank] * nHistoriesPerParticle],
              nPPerRank[world_rank] * nHistoriesPerParticle, MPI_FLOAT,
              &weightHistoryGather[0], phpn, displ, MPI_FLOAT, 0,
              MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  }

  if( spectroscopy > 0 )
  {
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&net_Bins[0], &net_BinsTotal[0], nSpec, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  }
  
  if( surface_model > 0 || flux_ea > 0 )
  {
  // MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "Starting surface reduce " << std::endl;
  // for(int i=0;i<nSurfaces;i++) std::cout <<
  // surfaces->grossDeposition[i]<<std::endl;
  MPI_Reduce(&surfaces->grossDeposition[0], &grossDeposition[0], nSurfaces,
             MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->grossErosion[0], &grossErosion[0], nSurfaces, MPI_FLOAT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->sumWeightStrike[0], &sumWeightStrike[0], nSurfaces,
             MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->aveSputtYld[0], &aveSputtYld[0], nSurfaces, MPI_FLOAT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->sputtYldCount[0], &sputtYldCount[0], nSurfaces, MPI_INT,
             MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->sumParticlesStrike[0], &sumParticlesStrike[0],
             nSurfaces, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->energyDistribution[0], &energyDistribution[0],
             nSurfaces * nEdist * nAdist, MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->sputtDistribution[0], &sputtDistribution[0],
             nSurfaces * nEdist * nAdist, MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&surfaces->reflDistribution[0], &reflDistribution[0],
             nSurfaces * nEdist * nAdist, MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout << "Finished surface reduce " << std::endl;
  }
#endif
  if (world_rank == 0) {
    auto MPIfinish_clock = gitr_time::now();
    std::chrono::duration<gitr_precision> fsmpi = MPIfinish_clock - finish_clock;
    printf("Time taken for mpi reduction          is %6.3f (secs) \n",
           fsmpi.count());
  }
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    int totalHitWall = 0;
    for (int i = 0; i < nP; i++) {
      if (particleArray->hitWall[i] > 0.0)
        totalHitWall++;
    }
    if( use_3d_geom > 0 )
    {
    gitr_precision meanTransitTime0 = 0.0;
    meanTransitTime0 = meanTransitTime0 / nP;
    int max_boundary = 0;
    gitr_precision max_impacts = 0.0;
    int max_boundary1 = 0;
    gitr_precision max_impacts1 = 0.0;
    gitr_precision *impacts = new gitr_precision[nLines];
    gitr_precision *xOut = new gitr_precision[nP];
    gitr_precision *redeposit = new gitr_precision[nLines];
    gitr_precision *startingParticles = new gitr_precision[nLines];
    gitr_precision *surfZ = new gitr_precision[nLines];

    for (int i = 0; i < nLines; i++) {
      impacts[i] = boundaries[i].impacts;
      redeposit[i] = boundaries[i].redeposit;
      startingParticles[i] = boundaries[i].startingParticles;
      if (boundaries[i].impacts > max_impacts) {
        max_impacts = boundaries[i].impacts;
        max_boundary = i;
      }
      surfZ[i] = boundaries[i].Z;
    }

    for (int i = 0; i < nP; i++) {
      xOut[i] = particleArray->x[i];
    }

    }
    else
    {
//#else
  gitr_precision *impacts = new gitr_precision[nLines];
  gitr_precision *startingParticles = new gitr_precision[nLines];
  gitr_precision *surfZ = new gitr_precision[nLines];
  // float* impactEnergy = new float[nLines*1000];
  for (int i = 0; i < nLines; i++) {
    impacts[i] = boundaries[i].impacts;
    startingParticles[i] = boundaries[i].startingParticles;
    surfZ[i] = boundaries[i].Z;
  }
//#endif
    }
    // add initial particle erosion to surface counting
    int closestBoundaryIndex = 0;
    int surfIndex = 0;
    gitr_precision minDistance = 0.0;
    gitr_precision thisE[3] = {0.0};
    for (int j = 0; j < nP; j++) {
      minDistance =
          getE(px[j], py[j], pz[j], thisE, boundaries.data(), nLines,
               nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
               n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
               &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
               &closeGeom_sheath.front(), closestBoundaryIndex, biased_surface,
               use_3d_geom, geom_hash_sheath, cylsymm, sheath_model_type);
      
      if (boundaries[closestBoundaryIndex].Z > 0.0) {
        surfIndex = boundaries[closestBoundaryIndex].surfaceNumber;
        grossErosion[surfIndex] = grossErosion[surfIndex] + 1.0;
      }
    }
    
    ofstream outfile2;
    outfile2.open("output/positions.m");
    for (int i = 1; i < nP + 1; i++) {
      outfile2 << "Pos( " << i << ",:) = [ ";
      outfile2 << particleArray->x[i - 1] << " " << particleArray->y[i - 1]
               << " " << particleArray->z[i - 1] << " ];" << std::endl;
    }
    outfile2.close();

    // Write netCDF output for positions
    netCDF::NcFile ncFile0("output/positions.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nP0 = ncFile0.addDim("nP", nP);
    vector<netCDF::NcDim> dims0;
    dims0.push_back(nc_nP0);

    netCDF::NcVar nc_x0 = ncFile0.addVar("x", netcdf_precision, dims0);
    netCDF::NcVar nc_y0 = ncFile0.addVar("y", netcdf_precision, dims0);
    netCDF::NcVar nc_z0 = ncFile0.addVar("z", netcdf_precision, dims0);
    netCDF::NcVar nc_vx0 = ncFile0.addVar("vx", netcdf_precision, dims0);
    netCDF::NcVar nc_vy0 = ncFile0.addVar("vy", netcdf_precision, dims0);
    netCDF::NcVar nc_vz0 = ncFile0.addVar("vz", netcdf_precision, dims0);
    netCDF::NcVar nc_trans0 = ncFile0.addVar("transitTime", netcdf_precision, dims0);
    netCDF::NcVar nc_impact0 = ncFile0.addVar("hitWall", netcdf_precision, dims0);
    netCDF::NcVar nc_surfHit0 = ncFile0.addVar("surfaceHit", netCDF::ncInt, dims0);
    netCDF::NcVar nc_weight0 = ncFile0.addVar("weight", netcdf_precision, dims0);
    netCDF::NcVar nc_charge0 = ncFile0.addVar("charge", netcdf_precision, dims0);
    netCDF::NcVar nc_leak0 = ncFile0.addVar("hasLeaked", netCDF::ncInt, dims0);
    netCDF::NcVar nc_dist0 = ncFile0.addVar("distTraveled", netcdf_precision, dims0);
    netCDF::NcVar nc_time0 = ncFile0.addVar("time", netcdf_precision, dims0);
    netCDF::NcVar nc_dt0 = ncFile0.addVar("dt", netcdf_precision, dims0);
    // add mass and Z
    netCDF::NcVar nc_mass0 = ncFile0.addVar("amu", netcdf_precision, dims0);
    netCDF::NcVar nc_Z0 = ncFile0.addVar("Z", netcdf_precision, dims0);
    // species type
    netCDF::NcVar nc_species0 = ncFile0.addVar("species", netCDF::ncInt, dims0);
#if USE_MPI > 0
    nc_x0.putVar(&xGather[0]);
    nc_y0.putVar(&yGather[0]);
    nc_z0.putVar(&zGather[0]);
    nc_vx0.putVar(&vxGather[0]);
    nc_vy0.putVar(&vyGather[0]);
    nc_vz0.putVar(&vzGather[0]);
    nc_trans0.putVar(&particleArray->transitTime[0]);
    nc_impact0.putVar(&hitWallGather[0]);
    nc_surfHit0.putVar(&surfaceHitGather[0]);
    nc_weight0.putVar(&weightGather[0]);
    nc_charge0.putVar(&chargeGather[0]);
    nc_leak0.putVar(&hasLeakedGather[0]);
#else
  std::cout << "not using mpi output" << std::endl;
  nc_x0.putVar(&particleArray->xprevious[0]);
  nc_y0.putVar(&particleArray->yprevious[0]);
  nc_z0.putVar(&particleArray->zprevious[0]);
  nc_vx0.putVar(&particleArray->vx[0]);
  nc_vy0.putVar(&particleArray->vy[0]);
  nc_vz0.putVar(&particleArray->vz[0]);
  nc_trans0.putVar(&particleArray->transitTime[0]);
  nc_impact0.putVar(&particleArray->hitWall[0]);
  nc_surfHit0.putVar(&particleArray->surfaceHit[0]);
  nc_weight0.putVar(&particleArray->weight[0]);
  nc_charge0.putVar(&particleArray->charge[0]);
  nc_leak0.putVar(&particleArray->hasLeaked[0]);
  nc_dist0.putVar(&particleArray->distTraveled[0]);
  nc_time0.putVar(&particleArray->time[0]);
  nc_dt0.putVar(&particleArray->dt[0]);
  nc_mass0.putVar(&particleArray->amu[0]);
  nc_Z0.putVar(&particleArray->Z[0]);
  nc_species0.putVar(&particleArray->species[0]);
#endif
    ncFile0.close();
  if( surface_model > 0 || flux_ea > 0 )
  {
//// FIXME -- dump surface file --> fix flattening arrays

#if USE_MPI > 0
    std::vector<int> surfaceNumbers(nSurfaces, 0);
    int srf = 0;
    for (int i = 0; i < nLines; i++) {
      if (boundaries[i].surface) {
        surfaceNumbers[srf] = i;

        srf = srf + 1;
      }
    }
    netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
    vector<netCDF::NcDim> dims1;
    dims1.push_back(nc_nLines);

    // get species names

    vector<netCDF::NcDim> dimsSurfE;
    dimsSurfE.push_back(nc_nLines);
    netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
    netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
    dimsSurfE.push_back(nc_nAngles);
    dimsSurfE.push_back(nc_nEnergies);
    netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike", netcdf_precision, nc_nLines);
    nc_grossDep.putVar(&grossDeposition[0]);
    nc_surfNum.putVar(&surfaceNumbers[0]);
    nc_grossEro.putVar(&grossErosion[0]);
    nc_aveSpyl.putVar(&aveSputtYld[0]);
    nc_spylCounts.putVar(&sputtYldCount[0]);
    nc_sumParticlesStrike.putVar(&sumParticlesStrike[0]);
    nc_sumWeightStrike.putVar(&sumWeightStrike[0]);
    netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
    netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
    netCDF::NcVar nc_surfSputtDist =  ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);
    nc_surfEDist.putVar(&energyDistribution[0]);
    nc_surfReflDist.putVar(&reflDistribution[0]);
    nc_surfSputtDist.putVar(&sputtDistribution[0]);
    ncFile1.close();
#else
    int nSpecies = 0;
    std::vector<int> surfaceNumbers(nSurfaces, 0);
    int srf = 0;
    for (int i = 0; i < nLines; i++) {
      if (boundaries[i].surface) {
        surfaceNumbers[srf] = i;

        surfaces->grossErosion[srf] = surfaces->grossErosion[srf] + grossErosion[srf];
        srf = srf + 1;
      }
    }  
    netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
    vector<netCDF::NcDim> dims1;
    dims1.push_back(nc_nLines);

    vector<netCDF::NcDim> dimsSurfE;
    dimsSurfE.push_back(nc_nLines);
    netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
    netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
    dimsSurfE.push_back(nc_nEnergies);
    dimsSurfE.push_back(nc_nAngles);
    netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, nc_nLines);
    netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, nc_nLines);
    netCDF::NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike", netcdf_precision, nc_nLines);
    nc_grossDep.putVar(&surfaces->grossDeposition[0]);
    nc_surfNum.putVar(&surfaceNumbers[0]);
    nc_grossEro.putVar(&surfaces->grossErosion[0]);
    nc_aveSpyl.putVar(&surfaces->aveSputtYld[0]);
    nc_spylCounts.putVar(&surfaces->sputtYldCount[0]);
    nc_sumParticlesStrike.putVar(&surfaces->sumParticlesStrike[0]);
    nc_sumWeightStrike.putVar(&surfaces->sumWeightStrike[0]);
    netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
    netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
    netCDF::NcVar nc_surfSputtDist = ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);
    nc_surfEDist.putVar(&surfaces->energyDistribution[0]);
    nc_surfReflDist.putVar(&surfaces->reflDistribution[0]);
    nc_surfSputtDist.putVar(&surfaces->sputtDistribution[0]);
    ncFile1.close();
  }
#endif
  if( particle_tracks > 0 )
  {
    // Write netCDF output for histories
    netCDF::NcFile ncFile_hist("output/history.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nT = ncFile_hist.addDim("nT", nHistoriesPerParticle);
    netCDF::NcDim nc_nP = ncFile_hist.addDim("nP", nP);
    vector<netCDF::NcDim> dims_hist;
    dims_hist.push_back(nc_nP);
    dims_hist.push_back(nc_nT);
    netCDF::NcVar nc_x = ncFile_hist.addVar("x", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_y = ncFile_hist.addVar("y", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_z = ncFile_hist.addVar("z", netCDF::ncDouble, dims_hist);

    netCDF::NcVar nc_v = ncFile_hist.addVar("v", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vx = ncFile_hist.addVar("vx", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vy = ncFile_hist.addVar("vy", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vz = ncFile_hist.addVar("vz", netCDF::ncDouble, dims_hist);

    netCDF::NcVar nc_charge = ncFile_hist.addVar("charge", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_Z = ncFile_hist.addVar("Z", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_weight = ncFile_hist.addVar("weight", netCDF::ncDouble, dims_hist);

#if USE_MPI > 0
    nc_x.putVar(&positionHistoryXgather[0]);
    nc_y.putVar(&positionHistoryYgather[0]);
    nc_z.putVar(&positionHistoryZgather[0]);

    nc_v.putVar(&velocityHistorygather[0]);
    nc_vx.putVar(&velocityHistoryXgather[0]);
    nc_vy.putVar(&velocityHistoryYgather[0]);
    nc_vz.putVar(&velocityHistoryZgather[0]);

    nc_charge.putVar(&chargeHistoryGather[0]);
    nc_Z.putVar(&ZHistoryGather[0]);
    nc_weight.putVar(&weightHistoryGather[0]);
#else
    nc_x.putVar(&positionHistoryX[0]);
    nc_y.putVar(&positionHistoryY[0]);
    nc_z.putVar(&positionHistoryZ[0]);

    nc_vx.putVar(&velocityHistoryX[0]);
    nc_vy.putVar(&velocityHistoryY[0]);
    nc_vz.putVar(&velocityHistoryZ[0]);

    nc_charge.putVar(&chargeHistory[0]);
    nc_Z.putVar(&ZHistory[0]);
#endif
    ncFile_hist.close();
  }
    if( spectroscopy > 0 )
    {
    // Write netCDF output for density data
    netCDF::NcFile ncFile("output/spec.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nBins = ncFile.addDim("nBins", nBins + 1);
    netCDF::NcDim nc_nR = ncFile.addDim("nR", net_nX);
    netCDF::NcDim nc_nY;
    if( spectroscopy > 2 )
    {
    nc_nY = ncFile.addDim("nY", net_nY);
    }

    netCDF::NcDim nc_nZ = ncFile.addDim("nZ", net_nZ);
    vector<netCDF::NcDim> dims;
    dims.push_back(nc_nBins);
    dims.push_back(nc_nZ);

    if( spectroscopy > 2 )
    {
    dims.push_back(nc_nY);
    }

    dims.push_back(nc_nR);

    netCDF::NcVar nc_n = ncFile.addVar("n", netcdf_precision, dims);
    netCDF::NcVar nc_gridR = ncFile.addVar("gridR", netcdf_precision, nc_nR);
    netCDF::NcVar nc_gridZ = ncFile.addVar("gridZ", netcdf_precision, nc_nZ);
    nc_gridR.putVar(&gridX_bins[0]);
    nc_gridZ.putVar(&gridZ_bins[0]);
    if( spectroscopy > 2 )
    {
    netCDF::NcVar nc_gridY = ncFile.addVar("gridY", netcdf_precision, nc_nY);
    nc_gridY.putVar(&gridY_bins[0]);
    }

#if USE_MPI > 0
    nc_n.putVar(&net_BinsTotal[0]);
#else
    nc_n.putVar(&net_Bins[0]);
#endif
    ncFile.close();
    }
#ifdef __CUDACC__
    cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
#endif
#if USE_MPI > 0
  }
#endif
#ifdef __CUDACC__
  cudaError_t err = cudaDeviceReset();
// cudaProfilerStop();
#endif
  if (world_rank == 0) {
    auto gitr_finish_clock = gitr_time::now();
    std::chrono::duration<gitr_precision> fstotal = gitr_finish_clock - gitr_start_clock;
    printf("Total runtime for GITR is %6.3f (secs) \n", fstotal.count());
  }
#if USE_MPI > 0
  // Finalize the MPI environment.
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return 0;
}

