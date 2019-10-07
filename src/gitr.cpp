#include "Boundary.h"
#include "Fields.h"
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
#include "interpRateCoeff.hpp"
#include "interpolate.h"
#include "ionize.h"
#include <cmath>
//#include "ncFile.h"
#include "ompPrint.h"
#include "recombine.h"
#include "spectroscopy.h"
#include "surfaceModel.h"
#include "testRoutineCuda.h"
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

#ifdef __CUDACC__
#include <curand.h>
#include <curand_kernel.h>
#else
#endif

#if USE_FILESYSTEM > 0
#include <experimental/filesystem>
#endif

#if USE_MPI
#include <mpi.h>
#endif

#include <omp.h>

#include "sortParticles.h"
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

int main(int argc, char **argv, char **envp) {
  typedef std::chrono::high_resolution_clock gitr_time;
  auto gitr_start_clock = gitr_time::now();

  //Set default processes per node to 1
  int ppn = 1;

  //Set default input file string
  std::string inputFile = "gitrInput.cfg";

#if USE_MPI > 0
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
#endif
 
  // read comand line arguments for specifying number of ppn (or gpus per node)
  // and specify input file if different than default
  // -nGPUPerNode and -i respectively
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
  
  // Prepare config files for import
  libconfig::Config cfg, cfg_geom;
  std::string input_path = "input/";
  
  if (world_rank == 0) {
    // Parse and read input file
    std::cout << "Open configuration file " << input_path << inputFile
              << std::endl;
    importLibConfig(cfg, input_path + inputFile);

    // Parse and read geometry file
    std::string geomFile;
    getVariable(cfg, "geometry.fileString", geomFile);
    std::cout << "Open geometry file " << input_path + geomFile << std::endl;
    importLibConfig(cfg_geom, input_path + geomFile);

    std::cout << "Successfully staged input and geometry file " << std::endl;

// check binary compatibility with input file
#if CHECK_COMPATIBILITY > 0
    checkFlags(cfg);
#endif
  }

// show memory usage of GPU
#if USE_FILESYSTEM > 0
#if __CUDACC__
  namespace fsn = std::experimental::filesystem;
#else
  namespace fsn = std::experimental::filesystem;
#endif

print_gpu_memory_usage(world_rank);

  fsn::path output_folder = "output";
  // Output

  //boost::filesystem::path dir(output_folder);
  if (!(fsn::exists(output_folder))) {
    std::cout << "Doesn't Exist in main" << std::endl;
    if (fsn::create_directory(output_folder)) {
      std::cout << " Successfully Created " << std::endl;
    }
  }
#endif
  // Background species info
  float background_Z = 0.0, background_amu = 0.0;
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
  typedef std::chrono::duration<float> fsec0nc;
  fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
  //printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
  
  
  int nR_Bfield = 1, nY_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  std::string bfieldFile;
  if (world_rank == 0) {
    importVectorFieldNs(cfg, input_path, BFIELD_INTERP, bfieldCfg, nR_Bfield,
                        nY_Bfield, nZ_Bfield, bfieldFile);
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Bfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  sim::Array<float> bfieldGridr(nR_Bfield), bfieldGridy(nY_Bfield),
      bfieldGridz(nZ_Bfield);
  n_Bfield = nR_Bfield * nY_Bfield * nZ_Bfield;
  sim::Array<float> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  if (world_rank == 0) {
    importVectorField(cfg, input_path, BFIELD_INTERP, bfieldCfg, nR_Bfield,
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
  float Btest[3] = {0.0f};
  interp2dVector(&Btest[0], 5.5, 0.0, -4.0, nR_Bfield, nZ_Bfield,
                 bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(),
                 by.data());
  //std::cout << "node " << world_rank << "Bfield at 5.5 -4 " << Btest[0] << " "
  //          << Btest[1] << " " << Btest[2] << std::endl;
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
    nSurfaces = importGeometry(cfg_geom, boundaries);
    std::cout << "Starting Boundary Init... nSurfaces " << nSurfaces
              << std::endl;
  }
#if USE_MPI > 0
  MPI_Bcast(&nSurfaces, 1, MPI_INT, 0, MPI_COMM_WORLD);
#if USE3DTETGEOM > 0
  const int nBoundaryMembers = 39;
  int nIntMembers = 5;
#else
  const int nBoundaryMembers = 37;
  int nIntMembers = 5;
#endif
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

  float biasPotential = 0.0;
#if BIASED_SURFACE > 0
  if (world_rank == 0) {
    getVariable(cfg, "backgroundPlasmaProfiles.biasPotential", biasPotential);
  }
#if USE_MPI > 0
  MPI_Bcast(&biasPotential, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  // create Surface data structures
  int nEdist = 1;
  float E0dist = 0.0;
  float Edist = 0.0;
  int nAdist = 1;
  float A0dist = 0.0;
  float Adist = 0.0;
#if FLUX_EA > 0
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
#endif
  auto surfaces = new Surfaces(nSurfaces, nEdist, nAdist);
  surfaces->setSurface(nEdist, E0dist, Edist, nAdist, A0dist, Adist);

  //#if USE_MPI > 0
  // Arrays used for reduction at end of sim
  sim::Array<float> grossDeposition(nSurfaces, 0.0);
  sim::Array<float> grossErosion(nSurfaces, 0.0);
  sim::Array<float> sumWeightStrike(nSurfaces, 0.0);
  sim::Array<float> energyDistribution(nSurfaces * nEdist * nAdist, 0.0);
  sim::Array<float> reflDistribution(nSurfaces * nEdist * nAdist, 0.0);
  sim::Array<float> sputtDistribution(nSurfaces * nEdist * nAdist, 0.0);
  sim::Array<float> aveSputtYld(nSurfaces, 0.0);
  sim::Array<int> sputtYldCount(nSurfaces, 0);
  sim::Array<int> sumParticlesStrike(nSurfaces, 0);
  //#endif

  int nHashes = 1;
  int nR_closeGeomTotal = 1;
  int nY_closeGeomTotal = 1;
  int nZ_closeGeomTotal = 1;
  int nHashPointsTotal = 1;
  int nGeomHash = 1;
  std::string geomHashCfg = "geometry_hash.";
#if GEOM_HASH == 1
  nR_closeGeomTotal = 0;
  nY_closeGeomTotal = 0;
  nZ_closeGeomTotal = 0;
  nHashPointsTotal = 0;
  nGeomHash = 0;
  if (world_rank == 0) {
    getVariable(cfg, geomHashCfg + "nHashes", nHashes);
  }
#if USE_MPI > 0
  MPI_Bcast(&nHashes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  sim::Array<int> nR_closeGeom(nHashes, 0);
  sim::Array<int> nY_closeGeom(nHashes, 0);
  sim::Array<int> nZ_closeGeom(nHashes, 0);
  sim::Array<int> nHashPoints(nHashes, 0);
  sim::Array<int> n_closeGeomElements(nHashes, 0);

#if GEOM_HASH == 1
  if (world_rank == 0) {
    importHashNs(cfg, input_path, nHashes, "geometry_hash", nR_closeGeom.data(),
                 nY_closeGeom.data(), nZ_closeGeom.data(),
                 n_closeGeomElements.data(), nR_closeGeomTotal,
                 nY_closeGeomTotal, nZ_closeGeomTotal, nHashPoints.data(),
                 nHashPointsTotal, nGeomHash);
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
    }
     for(int j=0;j<nHashes;j++)
    {
      nGeomHash = nGeomHash +
      nR_closeGeom[j]*nZ_closeGeom[j]*n_closeGeomElements[j];
      nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[j];
      nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[j];
    }
    #if USE3DTETGEOM > 0
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
    #endif
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
#endif

#if GEOM_HASH > 1
  if (world_rank == 0) {
    getVariable(cfg, geomHashCfg + "nHashes", nHashes);
  }
#if USE_MPI > 0
  MPI_Bcast(&nHashes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  std::vector<std::string> hashFile;
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
#if USE3DTETGEOM > 0
      nY_closeGeom[i] = getDimFromFile(cfg, input_path + hashFile[i],
                                       geomHashCfg, "gridNyString");
      nGeomHash = nGeomHash -
                  nR_closeGeom[i] * nZ_closeGeom[i] * n_closeGeomElements[i] +
                  nY_closeGeom[i] * nR_closeGeom[i] * nZ_closeGeom[i] *
                      n_closeGeomElements[i];
      nY_closeGeomTotal = nY_closeGeomTotal + nY_closeGeom[i];
#endif
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
#endif

  std::cout << "allocating closGeomGrids " << nR_closeGeomTotal << " "
            << nY_closeGeomTotal << " " << nZ_closeGeomTotal << " " << nGeomHash
            << std::endl;
  sim::Array<float> closeGeomGridr(nR_closeGeomTotal),
      closeGeomGridy(nY_closeGeomTotal), closeGeomGridz(nZ_closeGeomTotal);
  sim::Array<int> closeGeom(nGeomHash, 0);
  std::cout << "allocating closGeomGrids finished" << std::endl;

#if GEOM_HASH == 1
  sim::Array<float> hashX0(nHashes, 0.0), hashX1(nHashes, 0.0),
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
#if USE3DTETGEOM > 0
        hashY0[i] = geomHash["hashY0"][i];
        hashY1[i] = geomHash["hashY1"][i];
#endif
      }
    } else {
      getVariable(cfg, geomHashCfg + "hashX0", hashX0[0]);
      getVariable(cfg, geomHashCfg + "hashX1", hashX1[0]);
      getVariable(cfg, geomHashCfg + "hashZ0", hashZ0[0]);
      getVariable(cfg, geomHashCfg + "hashZ1", hashZ1[0]);
#if USE3DTETGEOM > 0
      getVariable(cfg, geomHashCfg + "hashY0", hashY0[0]);
      getVariable(cfg, geomHashCfg + "hashY1", hashY1[0]);
#endif
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
    // std::cout << "gridX "<< closeGeomGridr[i] << std::endl;
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
    // std::cout << "gridY "<< closeGeomGridy[j] << std::endl;
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
    // std::cout << "gridz "<< closeGeomGridz[k] << std::endl;
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
  sim::Array<float> minDist1(Maxn_closeGeomElements, 1e6);
  std::cout << "Generating geometry hash" << sizeof(int) << " bytes per int, "
            << nGeomHash << " for the entire hash " << std::endl;

#if USE_CUDA > 0
  // cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

  // if(cudaSuccess != cuda_status )
  //{

  //  printf("Error: cudaMemGetInfo fails, %s \n",
  //  cudaGetErrorString(cuda_status) ); exit(1);
  //}

  // free_db = (double)free_byte ;
  // total_db = (double)total_byte ;
  // used_db = total_db - free_db ;

  // printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
  //  used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
#endif
  std::cout << "starting geomhash" << std::endl;
  typedef std::chrono::high_resolution_clock Time0;
  typedef std::chrono::duration<float> fsec0;
  auto start_clock0 = Time0::now();
  hashGeom geo1(nLines, nHashes, boundaries.data(), closeGeomGridr.data(),
                closeGeomGridy.data(), closeGeomGridz.data(),
                n_closeGeomElements.data(), closeGeom.data(),
                nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data());
  thrust::for_each(thrust::device,
                   lines0 + world_rank * nHashMeshPointsPerProcess,
                   lines0 + world_rank * nHashMeshPointsPerProcess +
                       hashMeshIncrements[world_rank] - 1,
                   geo1);
// for(int i=0;i<nR_closeGeom*nY_closeGeom*nZ_closeGeom;i++)
//{
// geo1(i);
//}
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
#if USE3DTETGEOM
        for (int j = 0; j < nY_closeGeom[ii]; j++) {
#endif
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
#if USE3DTETGEOM
        }
#endif
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
      NcFile ncFile_hash("output/geomHash" + std::to_string(i) + ".nc",
                         NcFile::replace);
  std::cout << "opened GeomHash file " << std::endl;
      NcDim hashNR = ncFile_hash.addDim("nR", nR_closeGeom[i]);
  std::cout << "added first dimension nR " << std::endl;
#if USE3DTETGEOM > 0
      NcDim hashNY = ncFile_hash.addDim("nY", nY_closeGeom[i]);
#endif
      NcDim hashNZ = ncFile_hash.addDim("nZ", nZ_closeGeom[i]);
      NcDim hashN = ncFile_hash.addDim("n", n_closeGeomElements[i]);
      vector<NcDim> geomHashDim;
      geomHashDim.push_back(hashNR);
#if USE3DTETGEOM > 0
      geomHashDim.push_back(hashNY);
#endif
      geomHashDim.push_back(hashNZ);
      geomHashDim.push_back(hashN);
  std::cout << "added all dimensions " << std::endl;
      NcVar hash_gridR = ncFile_hash.addVar("gridR", ncFloat, hashNR);
  std::cout << "added first variable " << std::endl;
#if USE3DTETGEOM > 0
      NcVar hash_gridY = ncFile_hash.addVar("gridY", ncFloat, hashNY);
#endif
      NcVar hash_gridZ = ncFile_hash.addVar("gridZ", ncFloat, hashNZ);
  std::cout << "added gridZ variable " << std::endl;
      NcVar hash = ncFile_hash.addVar("hash", ncInt, geomHashDim);
  std::cout << "added hash variable " << std::endl;
      int ncIndex = 0;
      if (i > 0)
        ncIndex = nR_closeGeom[i - 1];
      hash_gridR.putVar(&closeGeomGridr[ncIndex]);
  std::cout << "put gridr variable " << std::endl;
#if USE3DTETGEOM > 0
      if (i > 0)
        ncIndex = nY_closeGeom[i - 1];
      hash_gridY.putVar(&closeGeomGridy[ncIndex]);
  std::cout << "put gridy variable " << std::endl;
#endif

      if (i > 0)
        ncIndex = nZ_closeGeom[i - 1];
      hash_gridZ.putVar(&closeGeomGridz[ncIndex]);
  std::cout << "put gridz variable " << std::endl;
      if (i > 0)
        ncIndex = nR_closeGeom[i - 1] * nY_closeGeom[i - 1] *
                  nZ_closeGeom[i - 1] * n_closeGeomElements[i - 1];
  std::cout << "ncIndex "<<ncIndex << std::endl;
      hash.putVar(&closeGeom[ncIndex]);
  std::cout << "put hash variable " << std::endl;
      ncFile_hash.close();
  std::cout << "closed file " << std::endl;
    }
  }
  std::cout << "finished GeomHash file " << std::endl;
#elif GEOM_HASH > 1
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
#if USE3DTETGEOM > 0
      if (i > 0)
        dataIndex = nY_closeGeom[0];
      getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg, "gridYString",
                     closeGeomGridy[dataIndex]);
#endif
      if (i > 0)
        dataIndex = nR_closeGeom[0] * nY_closeGeom[0] * nZ_closeGeom[0] *
                    n_closeGeomElements[0];
      getVarFromFile(cfg, input_path + hashFile[i], geomHashCfg,
                     "closeGeomString", closeGeom[dataIndex]);
    }
  }
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
#endif

  int nR_closeGeom_sheath = 1;
  int nY_closeGeom_sheath = 1;
  int nZ_closeGeom_sheath = 1;
  int n_closeGeomElements_sheath = 1;
  int nGeomHash_sheath = 1;
  std::string geomHashSheathCfg = "geometry_sheath.";
#if GEOM_HASH_SHEATH == 1
  if (world_rank == 0) {
    getVariable(cfg, geomHashSheathCfg + "nR_closeGeom", nR_closeGeom_sheath);
    getVariable(cfg, geomHashSheathCfg + "nZ_closeGeom", nZ_closeGeom_sheath);
    getVariable(cfg, geomHashSheathCfg + "n_closeGeomElements",
                n_closeGeomElements_sheath);
    nGeomHash_sheath =
        nR_closeGeom_sheath * nZ_closeGeom_sheath * n_closeGeomElements_sheath;
#if USE3DTETGEOM > 0
    getVariable(cfg, geomHashSheathCfg + "nY_closeGeom", nY_closeGeom_sheath);
    nGeomHash_sheath = nY_closeGeom_sheath * nGeomHash_sheath;
#endif
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

#if GEOM_HASH_SHEATH > 1
  std::string hashFile_sheath;
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
#if USE3DTETGEOM > 0
    nY_closeGeom_sheath = getDimFromFile(cfg, input_path + hashFile_sheath,
                                         geomHashSheathCfg, "gridNyString");
    nGeomHash_sheath = nY_closeGeom_sheath * nGeomHash_sheath;
#endif
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_closeGeom_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_closeGeomElements_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nGeomHash_sheath, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  sim::Array<float> closeGeomGridr_sheath(nR_closeGeom_sheath),
      closeGeomGridy_sheath(nY_closeGeom_sheath),
      closeGeomGridz_sheath(nZ_closeGeom_sheath);
  sim::Array<int> closeGeom_sheath(nGeomHash_sheath);

#if GEOM_HASH_SHEATH == 1
  float hashX0_s, hashX1_s, hashY0_s, hashY1_s, hashZ0_s, hashZ1_s;
  if (world_rank == 0) {
    getVariable(cfg, geomHashSheathCfg + "hashX0", hashX0_s);
    getVariable(cfg, geomHashSheathCfg + "hashX1", hashX1_s);
    getVariable(cfg, geomHashSheathCfg + "hashZ0", hashZ0_s);
    getVariable(cfg, geomHashSheathCfg + "hashZ1", hashZ1_s);
#if USE3DTETGEOM > 0
    getVariable(cfg, geomHashSheathCfg + "hashY0", hashY0_s);
    getVariable(cfg, geomHashSheathCfg + "hashY1", hashY1_s);
#endif
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
  sim::Array<float> minDist1_s(nGeomHash_sheath, 1e6);
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
  typedef std::chrono::duration<float> fsec0_s;
  auto start_clock0_s = Time0_s::now();
  hashGeom_sheath geo_s(
      nLines, boundaries.data(), closeGeomGridr_sheath.data(),
      closeGeomGridy_sheath.data(), closeGeomGridz_sheath.data(),
      n_closeGeomElements_sheath, closeGeom_sheath.data(), nR_closeGeom_sheath,
      nY_closeGeom_sheath, nZ_closeGeom_sheath);
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
    NcFile ncFile_hash_sheath("output/geomHash_sheath.nc", NcFile::replace);
    NcDim hashNR_sheath = ncFile_hash_sheath.addDim("nR", nR_closeGeom_sheath);
    NcDim hashNY_sheath = ncFile_hash_sheath.addDim("nY", nY_closeGeom_sheath);
    NcDim hashNZ_sheath = ncFile_hash_sheath.addDim("nZ", nZ_closeGeom_sheath);
    NcDim hashN_sheath =
        ncFile_hash_sheath.addDim("n", n_closeGeomElements_sheath);
    vector<NcDim> geomHashDim_sheath;
    geomHashDim_sheath.push_back(hashNR_sheath);
    geomHashDim_sheath.push_back(hashNY_sheath);
    geomHashDim_sheath.push_back(hashNZ_sheath);
    geomHashDim_sheath.push_back(hashN_sheath);
    NcVar hash_gridR_sheath =
        ncFile_hash_sheath.addVar("gridR", ncFloat, hashNR_sheath);
    NcVar hash_gridY_sheath =
        ncFile_hash_sheath.addVar("gridY", ncFloat, hashNY_sheath);
    NcVar hash_gridZ_sheath =
        ncFile_hash_sheath.addVar("gridZ", ncFloat, hashNZ_sheath);
    NcVar hash_sheath =
        ncFile_hash_sheath.addVar("hash", ncInt, geomHashDim_sheath);
    hash_gridR_sheath.putVar(&closeGeomGridr_sheath[0]);
    hash_gridY_sheath.putVar(&closeGeomGridy_sheath[0]);
    hash_gridZ_sheath.putVar(&closeGeomGridz_sheath[0]);
    hash_sheath.putVar(&closeGeom_sheath[0]);
    ncFile_hash_sheath.close();
  }
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#elif GEOM_HASH_SHEATH > 1
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridRString", closeGeomGridr_sheath[0]);
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridZString", closeGeomGridz_sheath[0]);
#if USE3DTETGEOM > 0
    getVarFromFile(cfg, input_path + hashFile_sheath, geomHashSheathCfg,
                   "gridYString", closeGeomGridy_sheath[0]);
#endif
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
#endif

  int nR_Lc = 1;
  int nY_Lc = 1;
  int nZ_Lc = 1;
  int nTracers = 1;
  float r0_Lc, r1_Lc, y0_Lc, y1_Lc, z0_Lc, z1_Lc, dr;
  int nTraceSteps;
  std::string connLengthCfg = "connectionLength.";
  std::string lcFile;
  if (world_rank == 0) {
#if GENERATE_LC > 0
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
#else
#if LC_INTERP > 1
    nR_Lc =
        getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNrString");
    nY_Lc =
        getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNyString");
    nZ_Lc =
        getDimFromFile(cfg, input_path + lcFile, connLengthCfg, "gridNzString");
#endif
#endif
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

#if USE3DTETGEOM > 0
  nTracers = nR_Lc * nY_Lc * nZ_Lc;
#else
  nTracers = nR_Lc * nZ_Lc;
#endif

  sim::Array<float> Lc(nTracers), s(nTracers);
  sim::Array<float> gridRLc(nR_Lc), gridYLc(nY_Lc), gridZLc(nZ_Lc);
  sim::Array<int> noIntersectionNodes(nTracers);
#if GENERATE_LC > 0
  float lcBuffer = 0.0;
  // if( !boost::filesystem::exists( lcFile ) )
  // {
  //   std::cout << "No pre-existing connection length file found" << std::endl;
#if USE3DTETGEOM > 0
  float dy_Lc = (y1_Lc - y0_Lc) / (nY_Lc - 1);
  for (int j = 0; j < nY_Lc; j++) {
    gridYLc[j] = y0_Lc + j * dy_Lc;
  }
#endif
  float dr_Lc = (r1_Lc - r0_Lc) / (nR_Lc - 1);
  for (int i = 0; i < nR_Lc; i++) {
    gridRLc[i] = r0_Lc + i * dr_Lc;
  }

  float dz_Lc = (z1_Lc - z0_Lc) / (nZ_Lc - 1);
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
  auto forwardTracerParticles = new Particles(nTracers);
  auto backwardTracerParticles = new Particles(nTracers);
  int addIndex = 0;
  std::cout << "Initializing tracer particles" << std::endl;

  for (int i = 0; i < nR_Lc; i++) {
    for (int j = 0; j < nY_Lc; j++) {
      for (int k = 0; k < nZ_Lc; k++) {
#if USE3DTETGEOM > 0
        addIndex = i + j * nR_Lc + k * nR_Lc * nY_Lc;
#else
        addIndex = i + k * nR_Lc;
#endif
        forwardTracerParticles->setParticle(addIndex, gridRLc[i], gridYLc[j],
                                            gridZLc[k], 0.0, 0.0, 0.0, 0, 0.0,
                                            0.0);
        backwardTracerParticles->setParticle(addIndex, gridRLc[i], gridYLc[j],
                                             gridZLc[k], 0.0, 0.0, 0.0, 0, 0.0,
                                             0.0);
      }
    }
  }

  // dummy surfaces for Lcs calculation (geometry_check)
  auto dummy_surfaces = new Surfaces(1, 1, 1);
  dummy_surfaces->setSurface(1, 1, 1, 1, 1, 1);

  typedef std::chrono::high_resolution_clock Time_trace;
  typedef std::chrono::duration<float> fsec_trace;
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
                                      &bz.front(), &by.front()));

    thrust::for_each(thrust::device, lcBegin, lcEnd,
                     field_line_trace(-1.0, backwardTracerParticles, dr,
                                      boundaries.data(), nLines, nR_Lc, nZ_Lc,
                                      gridRLc.data(), gridZLc.data(), Lc.data(),
                                      nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                                      &bfieldGridz.front(), &br.front(),
                                      &bz.front(), &by.front()));

    thrust::for_each(
        thrust::device, lcBegin, lcEnd,
        geometry_check(forwardTracerParticles, nLines, &boundaries[0],
                       dummy_surfaces, dr, ii, nR_closeGeom.data(),
                       nY_closeGeom.data(), nZ_closeGeom.data(),
                       n_closeGeomElements.data(), &closeGeomGridr.front(),
                       &closeGeomGridy.front(), &closeGeomGridz.front(),
                       &closeGeom.front(), 0, 0.0, 0.0, 0, 0.0, 0.0));

    thrust::for_each(
        thrust::device, lcBegin, lcEnd,
        geometry_check(backwardTracerParticles, nLines, &boundaries[0],
                       dummy_surfaces, dr, ii, nR_closeGeom.data(),
                       nY_closeGeom.data(), nZ_closeGeom.data(),
                       n_closeGeomElements.data(), &closeGeomGridr.front(),
                       &closeGeomGridy.front(), &closeGeomGridz.front(),
                       &closeGeom.front(), 0, 0.0, 0.0, 0, 0.0, 0.0));
  }
  auto finish_clock_trace = Time_trace::now();
  fsec_trace fstrace = finish_clock_trace - start_clock_trace;
  printf("Time taken          is %6.3f (secs) \n", fstrace.count());
  printf("Time taken per step is %6.3f (secs) \n",
         fstrace.count() / (float)nTraceSteps);
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
  sim::Array<float> forwardHitWall(nTracers, 0.0),
      backwardHitWall(nTracers, 0.0), forwardTracerX(nTracers, 0.0),
      backwardTracerX(nTracers, 0.0);
  sim::Array<float> forwardTracerY(nTracers, 0.0),
      backwardTracerY(nTracers, 0.0);
  sim::Array<float> forwardTracerZ(nTracers, 0.0),
      backwardTracerZ(nTracers, 0.0);
  sim::Array<float> forwardDistanceTraveled(nTracers, 0.0),
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
    float forwardDist = 0.0;
    float backwardDist = 0.0;
    for (int i = 0; i < nR_Lc; i++) {
      for (int j = 0; j < nY_Lc; j++) {
        for (int k = 0; k < nZ_Lc; k++) {
          // std::cout << "hitwall of tracers " <<
          // forwardTracerParticles->hitWall[addIndex] << std::endl;
#if USE3DTETGEOM > 0
          addIndex = i + j * nR_Lc + k * nR_Lc * nY_Lc;
#else
        addIndex = i + k * nR_Lc;
#endif
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
          //       if(forwardTracerParticles->distanceTraveled[addIndex] >
          //               backwardTracerParticles->distanceTraveled[addIndex])
          // if(forwardTracerParticles->distanceTraveled[addIndex] >
          //        0.5*Lc[addIndex])
          //{
          //  s[addIndex] =
          //  -(0.5*Lc[addIndex]-backwardTracerParticles->distanceTraveled[addIndex]);
          //}
          // else
          //{
          //  s[addIndex] =
          //  (0.5*Lc[addIndex]-forwardTracerParticles->distanceTraveled[addIndex]);
          //}
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

    NcFile ncFileLC("LcS.nc", NcFile::replace);
    vector<NcDim> dims_lc;
    NcDim nc_nTracers = ncFileLC.addDim("nTracers", nTracers);
    NcDim nc_nRLc = ncFileLC.addDim("nR", nR_Lc);
    dims_lc.push_back(nc_nRLc);

#if USE3DTETGEOM
    NcDim nc_nYLc = ncFileLC.addDim("nY", nY_Lc);
    dims_lc.push_back(nc_nYLc);
#endif

    NcDim nc_nZLc = ncFileLC.addDim("nZ", nZ_Lc);
    dims_lc.push_back(nc_nZLc);

    NcVar nc_Lc = ncFileLC.addVar("Lc", ncFloat, dims_lc);
    NcVar nc_s = ncFileLC.addVar("s", ncFloat, dims_lc);
    NcVar nc_ftx = ncFileLC.addVar("fx", ncFloat, dims_lc);
    NcVar nc_fty = ncFileLC.addVar("fy", ncFloat, dims_lc);
    NcVar nc_ftz = ncFileLC.addVar("fz", ncFloat, dims_lc);
    NcVar nc_btx = ncFileLC.addVar("bx", ncFloat, dims_lc);
    NcVar nc_bty = ncFileLC.addVar("by", ncFloat, dims_lc);
    NcVar nc_btz = ncFileLC.addVar("bz", ncFloat, dims_lc);
    NcVar nc_nI = ncFileLC.addVar("noIntersection", ncFloat, dims_lc);
    NcVar nc_gridRLc = ncFileLC.addVar("gridR", ncFloat, nc_nRLc);
#if USE3DTETGEOM
    NcVar nc_gridYLc = ncFileLC.addVar("gridY", ncFloat, nc_nYLc);
#endif
    NcVar nc_gridZLc = ncFileLC.addVar("gridZ", ncFloat, nc_nZLc);

    nc_Lc.putVar(&Lc[0]);
    nc_s.putVar(&s[0]);
    nc_ftx.putVar(&forwardTracerX[0]);
    nc_fty.putVar(&forwardTracerY[0]);
    nc_ftz.putVar(&forwardTracerZ[0]);
    nc_btx.putVar(&backwardTracerX[0]);
    nc_bty.putVar(&backwardTracerY[0]);
    nc_btz.putVar(&backwardTracerZ[0]);
    nc_nI.putVar(&noIntersectionNodes[0]);
    nc_gridRLc.putVar(&gridRLc[0]);
#if USE3DTETGEOM
    nc_gridYLc.putVar(&gridYLc[0]);
#endif
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
#endif

#if LC_INTERP > 1
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
#endif

  // Background Plasma Temperature Initialization
  int nR_Temp = 1;
  int nY_Temp = 1;
  int nZ_Temp = 1;
  int n_Temp = 1;
  std::string tempCfg = "backgroundPlasmaProfiles.Temperature.";
#if TEMP_INTERP > 0
  std::string tempFile;
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, tempCfg + "fileString", tempFile);
    nR_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNrString");
#if TEMP_INTERP > 1
    nZ_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNzString");
#endif
#if TEMP_INTERP > 2
    nY_Temp =
        getDimFromFile(cfg, input_path + tempFile, tempCfg, "gridNyString");
#endif
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  sim::Array<float> TempGridr(nR_Temp), TempGridz(nZ_Temp), TempGridy(nY_Temp);
  n_Temp = nR_Temp * nY_Temp * nZ_Temp;
  sim::Array<float> ti(n_Temp), te(n_Temp);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
#if TEMP_INTERP == 0
    getVariable(cfg, tempCfg + "ti", ti[0]);
    getVariable(cfg, tempCfg + "te", te[0]);
#else
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridRString",
                 TempGridr[0]);
#if TEMP_INTERP > 1
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridZString",
                 TempGridz[0]);
#endif
#if TEMP_INTERP > 2
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "gridYString",
                 TempGridy[0]);
#endif
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "IonTempString", ti[0]);
  getVarFromFile(cfg, input_path + tempFile, tempCfg, "ElectronTempString",
                 te[0]);
#endif
#if USE_MPI > 0
  }

  MPI_Bcast(TempGridr.data(), nR_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(TempGridy.data(), nY_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(TempGridz.data(), nZ_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ti.data(), n_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(te.data(), n_Temp, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  float testVec = 0.0;
  testVec = interp2dCombined(0.0, 0.1, 0.0, nR_Temp, nZ_Temp, TempGridr.data(),
                             TempGridz.data(), ti.data());
  std::cout << "Finished Temperature import " << testVec << std::endl;

  // Background Plasma Density Initialization
  int nR_Dens = 1;
  int nY_Dens = 1;
  int nZ_Dens = 1;
  int n_Dens = 1;
  std::string densCfg = "backgroundPlasmaProfiles.Density.";
#if DENSITY_INTERP > 0
  std::string densFile;
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, densCfg + "fileString", densFile);
    nR_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNrString");
#if DENSITY_INTERP > 1
    nZ_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNzString");
#endif
#if DENSITY_INTERP > 2
    nY_Dens =
        getDimFromFile(cfg, input_path + densFile, densCfg, "gridNyString");
#endif
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_Dens, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  sim::Array<float> DensGridr(nR_Dens), DensGridz(nZ_Dens), DensGridy(nY_Dens);
  n_Dens = nR_Dens * nY_Dens * nZ_Dens;
  sim::Array<float> ni(n_Dens), ne(n_Dens);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
#if DENSITY_INTERP == 0
    getVariable(cfg, densCfg + "ni", ni[0]);
    getVariable(cfg, densCfg + "ne", ne[0]);
#else
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridRString",
                 DensGridr[0]);
#if DENSITY_INTERP > 1
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridZString",
                 DensGridz[0]);
#endif
#if DENSITY_INTERP > 2
  getVarFromFile(cfg, input_path + densFile, densCfg, "gridYString",
                 DensGridy[0]);
#endif
  getVarFromFile(cfg, input_path + densFile, densCfg, "IonDensityString",
                 ni[0]);
  getVarFromFile(cfg, input_path + densFile, densCfg, "ElectronDensityString",
                 ne[0]);
#endif
    std::cout << "Finished density import "
              << interp2dCombined(5.5, 0.0, -4.4, nR_Dens, nZ_Dens,
                                  &DensGridr.front(), &DensGridz.front(),
                                  &ne.front())
              << " "
              << interp2dCombined(0.0, 0.1, 0.0, nR_Dens, nZ_Dens,
                                  &DensGridr.front(), &DensGridz.front(),
                                  &ne.front())
              << std::endl;
// for(int i=0;i<100;i++)
//{
//    std::cout << i*0.001 << " " <<
//    interp2dCombined(0.001*i,0.0,0.0,nR_Dens,nZ_Dens,
//                                     &DensGridr.front(),&DensGridz.front(),&ne.front())
//                                     << std::endl;
//}
// std::cout << " z=0.1" << std::endl;
// for(int i=0; i<100;i++)
//{
//    std::cout << i*0.001 << " " <<
//    interp2dCombined(0.001*i,0.0,0.1,nR_Dens,nZ_Dens,
//                                     &DensGridr.front(),&DensGridz.front(),&ne.front())
//                                     << std::endl;
//}
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
#if FLOWV_INTERP == 1
  nR_flowV = nR_Lc;
  nY_flowV = nY_Lc;
  nZ_flowV = nZ_Lc;
#endif
#if FLOWV_INTERP > 1
  std::string flowVFile;
#if USE_MPI > 0
  if (world_rank == 0) {
    getVariable(cfg, flowVCfg + "fileString", flowVFile);
    nR_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNrString");
    nZ_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNzString");
#endif
#if FLOWV_INTERP > 2
    nY_flowV =
        getDimFromFile(cfg, input_path + flowVFile, flowVCfg, "gridNyString");
#endif
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_flowV, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

  sim::Array<float> flowVGridr(nR_flowV), flowVGridy(nY_flowV),
      flowVGridz(nZ_flowV);
  n_flowV = nR_flowV * nY_flowV * nZ_flowV;
  sim::Array<float> flowVr(n_flowV), flowVz(n_flowV), flowVt(n_flowV);

#if USE_MPI > 0
  if (world_rank == 0) {
#endif
#if FLOWV_INTERP == 0
    getVariable(cfg, flowVCfg + "flowVr", flowVr[0]);
    getVariable(cfg, flowVCfg + "flowVy", flowVt[0]);
    getVariable(cfg, flowVCfg + "flowVz", flowVz[0]);
#else
#if FLOWV_INTERP > 1
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
#endif
#if FLOWV_INTERP > 2
  getVarFromFile(cfg, input_path + flowVFile, flowVCfg, "gridYString",
                 flowVGridy[0]);
#endif
#endif
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
#if FLOWV_INTERP == 1
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
  int nFlowVs = nR_Lc * nZ_Lc;
#if LC_INTERP == 3
  for (int i = 0; i < nY_flowV; i++)
    flowVGridy[i] = gridYLc[i];
  nFlowVs = nR_Lc * nY_Lc * nZ_Lc;
#endif
  float thisY = 0.0;
  float cs0 = 0.0;
  float teLocal = 0.0;
  float tiLocal = 0.0;
  float BLocal[3] = {0.0, 0.0, 0.0};
  float Bnorm[3] = {0.0, 0.0, 0.0};
  float Bmag = 0.0;
  int index = 0;
  float cs = 0.0;
  float absS = 0.0;
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
                                   &TempGridz.front(), &te.front());
        tiLocal = interp2dCombined(flowVGridr[i], thisY, flowVGridz[j], nR_Temp,
                                   nZ_Temp, &TempGridr.front(),
                                   &TempGridz.front(), &ti.front());
        cs0 =
            std::sqrt((teLocal + tiLocal) * 1.602e-19 / (background_amu * 1.66e-27));
        interp2dVector(&BLocal[0], flowVGridr[i], thisY, flowVGridz[j],
                       nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                       bfieldGridz.data(), br.data(), bz.data(), by.data());
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
  sim::Array<float> flowVrSub(nFlowVs), flowVzSub(nFlowVs), flowVySub(nFlowVs);
  sim::Array<int> noIntersectionNearestMax(nFlowVs);
  float surroundingMinimumR = 0.0;
  float surroundingMinimumY = 0.0;
  float surroundingMinimumZ = 0.0;
  int iterIndex = 0;
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
  NcVar nc_flowVr = ncFileFlow.addVar("flowVr", ncFloat, dimsFlowV);
  NcVar nc_flowVt = ncFileFlow.addVar("flowVt", ncFloat, dimsFlowV);
  NcVar nc_flowVz = ncFileFlow.addVar("flowVz", ncFloat, dimsFlowV);
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
#endif

  // Background plasma temperature gradient field intitialization
  int nR_gradT = 1;
  int nY_gradT = 1;
  int nZ_gradT = 1;
  int n_gradT = 1;
  std::string gradTCfg = "backgroundPlasmaProfiles.gradT.";
  std::string gradTFile;
  if (world_rank == 0) {
#if GRADT_INTERP > 0
    getVariable(cfg, gradTCfg + "fileString", gradTFile);
    nR_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNrString");
#endif
#if GRADT_INTERP > 1
    nZ_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNzString");
#endif
#if GRADT_INTERP > 2
    nY_gradT =
        getDimFromFile(cfg, input_path + gradTFile, gradTCfg, "gridNyString");
#endif
  }
#if USE_MPI > 0
  MPI_Bcast(&nR_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_gradT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  n_gradT = nR_gradT * nY_gradT * nZ_gradT;
  sim::Array<float> gradTGridr(nR_gradT), gradTGridy(nY_gradT),
      gradTGridz(nZ_gradT);
  sim::Array<float> gradTeR(n_gradT), gradTeZ(n_gradT), gradTeY(n_gradT),
      gradTiR(n_gradT), gradTiZ(n_gradT), gradTiY(n_gradT);

  if (world_rank == 0) {
#if GRADT_INTERP == 0
    getVariable(cfg, gradTCfg + "gradTeR", gradTeR[0]);
    getVariable(cfg, gradTCfg + "gradTeY", gradTeY[0]);
    getVariable(cfg, gradTCfg + "gradTeZ", gradTeZ[0]);
    getVariable(cfg, gradTCfg + "gradTiR", gradTiR[0]);
    getVariable(cfg, gradTCfg + "gradTiY", gradTiY[0]);
    getVariable(cfg, gradTCfg + "gradTiZ", gradTiZ[0]);
#else
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridRString",
                   gradTGridr[0]);
#if GRADT_INTERP > 1
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridZString",
                   gradTGridz[0]);
#endif
#if GRADT_INTERP > 2
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gridYString",
                   gradTGridy[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTeYString",
                   gradTeY[0]);
    getVarFromFile(cfg, input_path + gradTFile, gradTCfg, "gradTiYString",
                   gradTiY[0]);
#endif

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
#endif
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
  float gradTi[3] = {0.0};
  interp2dVector(&gradTi[0], 1.45, 0.0, -1.2, nR_gradT, nZ_gradT,
                 gradTGridr.data(), gradTGridz.data(), gradTiR.data(),
                 gradTiZ.data(), gradTiY.data());
  std::cout << "thermal gradient interpolation gradTi " << gradTi[0] << " "
            << gradTi[1] << " " << gradTi[2] << " " << std::endl;

  // Initialization of ionization and recombination coefficients
  int nCS_Ionize = 1, nCS_Recombine = 1;
  const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
      *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
      *recombDensGrid, *recombTempGrid, *recombRCvarChar;
  std::string ionizeFile, recombFile;
  int nTemperaturesIonize = 1, nDensitiesIonize = 1;
  int i0, i1, i2, i3, i4;
  int nTemperaturesRecombine = 1, nDensitiesRecombine = 1;
#if USEIONIZATION > 0
  if (world_rank == 0) {
    if (cfg.lookupValue("impurityParticleSource.ionization.fileString",
                        ionizeFile) &&
        cfg.lookupValue("impurityParticleSource.ionization.nChargeStateString",
                        ionizeNcs) &&
        cfg.lookupValue("impurityParticleSource.ionization.DensGridString",
                        ionizeNDens) &&
        cfg.lookupValue("impurityParticleSource.ionization.TempGridString",
                        ionizeNTemp) &&
        cfg.lookupValue("impurityParticleSource.ionization.TempGridVarName",
                        ionizeTempGrid) &&
        cfg.lookupValue("impurityParticleSource.ionization.DensGridVarName",
                        ionizeDensGrid) &&
        cfg.lookupValue("impurityParticleSource.ionization.CoeffVarName",
                        ionizeRCvarChar)) {
      std::cout << "Ionization rate coefficient file: " << ionizeFile
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get ionization string info from input file "
          << std::endl;
    }
#endif
#if USERECOMBINATION > 0
    if (cfg.lookupValue("impurityParticleSource.recombination.fileString",
                        recombFile) &&
        cfg.lookupValue(
            "impurityParticleSource.recombination.nChargeStateString",
            recombNcs) &&
        cfg.lookupValue("impurityParticleSource.recombination.DensGridString",
                        recombNDens) &&
        cfg.lookupValue("impurityParticleSource.recombination.TempGridString",
                        recombNTemp) &&
        cfg.lookupValue("impurityParticleSource.recombination.TempGridVarName",
                        recombTempGrid) &&
        cfg.lookupValue("impurityParticleSource.recombination.DensGridVarName",
                        recombDensGrid) &&
        cfg.lookupValue("impurityParticleSource.recombination.CoeffVarName",
                        recombRCvarChar)) {
      std::cout << "Recombination rate coefficient file: " << recombFile
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get ionization string info from input file "
          << std::endl;
    }
#endif
#if USEIONIZATION > 0
    i0 = read_profileNs(input_path + ionizeFile, ionizeNcs, recombNcs,
                        nCS_Ionize, nCS_Recombine);

    i1 = read_profileNs(input_path + ionizeFile, ionizeNDens, ionizeNTemp,
                        nDensitiesIonize, nTemperaturesIonize);

    i3 = read_profileNs(input_path + recombFile, recombNDens, recombNTemp,
                        nDensitiesRecombine, nTemperaturesRecombine);
  }
#if USE_MPI > 0
  MPI_Bcast(&nCS_Ionize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nTemperaturesIonize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nDensitiesIonize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nCS_Recombine, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nTemperaturesRecombine, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nDensitiesRecombine, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  sim::Array<float> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
                                         nDensitiesIonize);
  sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
      gridDensity_Ionization(nDensitiesIonize);
  sim::Array<float> rateCoeff_Recombination(
      nCS_Recombine * nTemperaturesRecombine * nDensitiesRecombine);
  sim::Array<float> gridTemperature_Recombination(nTemperaturesRecombine),
      gridDensity_Recombination(nDensitiesRecombine);
  if (world_rank == 0) {
#if USEIONIZATION > 0
    i2 = read_profiles(
        input_path + ionizeFile, nTemperaturesIonize, nDensitiesIonize,
        ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
        gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
#endif
#if USERECOMBINATION > 0
    i4 = read_profiles(
        input_path + recombFile, nTemperaturesRecombine, nDensitiesRecombine,
        recombTempGrid, gridTemperature_Recombination, recombDensGrid,
        gridDensity_Recombination, recombRCvarChar, rateCoeff_Recombination);
#endif
  }
#if USE_MPI > 0
  MPI_Bcast(&rateCoeff_Ionization[0],
            nCS_Ionize * nTemperaturesIonize * nDensitiesIonize, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&gridTemperature_Ionization[0], nTemperaturesIonize, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&gridDensity_Ionization[0], nDensitiesIonize, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&rateCoeff_Recombination[0],
            nCS_Recombine * nTemperaturesRecombine * nDensitiesRecombine,
            MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gridTemperature_Recombination[0], nTemperaturesRecombine,
            MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gridDensity_Recombination[0], nDensitiesRecombine, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Applying background values at material boundaries
  std::for_each(boundaries.begin(), boundaries.end() - 1,
                boundary_init(background_Z, background_amu, nR_Dens, nZ_Dens,
                              DensGridr.data(), DensGridz.data(), ni.data(),
                              ne.data(), nR_Bfield, nZ_Bfield,
                              bfieldGridr.data(), bfieldGridz.data(), br.data(),
                              bz.data(), by.data(), nR_Temp, nZ_Temp,
                              TempGridr.data(), TempGridz.data(), ti.data(),
                              te.data(), biasPotential));

  std::cout << "Completed Boundary Init " << std::endl;

  // Efield
  int nR_PreSheathEfield = 1;
  int nY_PreSheathEfield = 1;
  int nZ_PreSheathEfield = 1;
  int nPSEs = 1;
  std::string PSECfg = "backgroundPlasmaProfiles.Efield.";
// sim::Array<float> preSheathEGridy(1);
#if USEPRESHEATHEFIELD > 0

  std::cout << "Using presheath Efield " << std::endl;
#if PRESHEATH_INTERP == 1
  nR_PreSheathEfield = nR_Lc;
  nY_PreSheathEfield = nY_Lc;
  nZ_PreSheathEfield = nZ_Lc;
#endif
#if PRESHEATH_INTERP > 1
  std::string efieldFile;
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, PSECfg + "fileString", efieldFile);
    nR_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNrString");
    nZ_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNzString");
#if PRESHEATH_INTERP > 2
    nY_PreSheathEfield =
        getDimFromFile(cfg, input_path + efieldFile, PSECfg, "gridNyString");
#endif
#if USE_MPI > 0
  }
  MPI_Bcast(&nR_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nY_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nZ_PreSheathEfield, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;
  sim::Array<float> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
  sim::Array<float> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
#if PRESHEATH_INTERP == 0
    getVariable(cfg, PSECfg + "Er", PSEr[0]);
    getVariable(cfg, PSECfg + "Et", PSEt[0]);
    getVariable(cfg, PSECfg + "Ez", PSEz[0]);
#elif PRESHEATH_INTERP > 1
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridRString",
                 preSheathEGridr[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridZString",
                 preSheathEGridz[0]);
#if PRESHEATH_INTERP > 2
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "gridYString",
                 preSheathEGridy[0]);
#endif

  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "radialComponentString",
                 PSEr[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg,
                 "toroidalComponentString", PSEt[0]);
  getVarFromFile(cfg, input_path + efieldFile, PSECfg, "axialComponentString",
                 PSEz[0]);
#endif
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

#if PRESHEATH_INTERP == 1

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
  float teLocal1 = 0.0;
  float BLocal1[3] = {0.0, 0.0, 0.0};
  float Bnorm1[3] = {0.0, 0.0, 0.0};
  float Bmag1 = 0.0;
  int index1 = 0;
  float absS1 = 0.0;
  float Epar = 0.0;
  for (int i = 0; i < nR_PreSheathEfield; i++) {
#if LC_INTERP == 3
    for (int k = 0; k < nY_PreSheathEfield; k++) {
      thisY = preSheathEGridy[k];
#endif
      for (int j = 0; j < nZ_PreSheathEfield; j++) {
        teLocal1 = interp2dCombined(preSheathEGridr[i], 0.0, preSheathEGridz[j],
                                    nR_Temp, nZ_Temp, &TempGridr.front(),
                                    &TempGridz.front(), &te.front());
        interp2dVector(&BLocal1[0], gridRLc[i], 0.0, gridZLc[j], nR_Bfield,
                       nZ_Bfield, bfieldGridr.data(), bfieldGridz.data(),
                       br.data(), bz.data(), by.data());
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
  sim::Array<float> PSErSub(nPSEs), PSEzSub(nPSEs), PSEySub(nPSEs);

  for (int i = 0; i < nR_Lc; i++) {
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
  NcVar nc_PSEr = ncFileLC.addVar("PSEr", ncFloat, nc_nTracers);
  NcVar nc_PSEt = ncFileLC.addVar("PSEt", ncFloat, nc_nTracers);
  NcVar nc_PSEz = ncFileLC.addVar("PSEz", ncFloat, nc_nTracers);
  nc_PSEr.putVar(&PSEr[0]);
  nc_PSEt.putVar(&PSEt[0]);
  nc_PSEz.putVar(&PSEz[0]);
#endif
#else
  nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;
  sim::Array<float> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
  sim::Array<float> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);

#endif
  std::string outnamePSEfieldR = "PSEfieldR.m";
  std::string outnamePSEfieldZ = "PSEfieldZ.m";
  std::string outnamePSEGridR = "PSEgridR.m";
  std::string outnamePSEGridZ = "PSEgridZ.m";
  // OUTPUT1d(profiles_folder,outnamePSEGridR, nR_PreSheathEfield,
  // &preSheathEGridr.front()); OUTPUT1d(profiles_folder,outnamePSEGridZ,
  // nZ_PreSheathEfield, &preSheathEGridz.front());
  //
  // OUTPUT3d(profiles_folder,outnamePSEfieldR,
  // nR_PreSheathEfield,nY_PreSheathEfield, nZ_PreSheathEfield, &PSEr.front());
  // OUTPUT3d(profiles_folder,outnamePSEfieldZ,
  // nR_PreSheathEfield,nY_PreSheathEfield, nZ_PreSheathEfield, &PSEz.front());
  std::cout << "Completed presheath Efield Init " << std::endl;
  sim::Array<float> Efieldr(nR_Bfield * nZ_Bfield),
      Efieldz(nR_Bfield * nZ_Bfield), Efieldt(nR_Bfield * nZ_Bfield),
      minDist(nR_Bfield * nZ_Bfield);

#if USESHEATHEFIELD > 0
#if EFIELD_INTERP == 1
  float thisE[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nR_Bfield; i++) {
    for (int j = 0; j < nZ_Bfield; j++) {
      minDist[(nR_Bfield - 1 - i) * nZ_Bfield + (nZ_Bfield - 1 - j)] =
          getE(bfieldGridr[i], 0.0, bfieldGridz[j], thisE, boundaries.data(),
               nLines, closestBoundaryIndex);
      Efieldr[i * nZ_Bfield + j] = thisE[0];
      Efieldz[i * nZ_Bfield + j] = thisE[2];
      Efieldt[i * nZ_Bfield + j] = thisE[1];
    }
  }

  int nR_closeGeom;
  int nZ_dtsEfield;

  int d1 = read_profileNs(
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNrString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNzString"),
      nR_dtsEfield, nZ_dtsEfield);

  sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
  sim::Array<float> dtsE(nR_dtsEfield * nZ_dtsEfield);

  int d2 = read_profile1d(
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridRString"),
      dtsEfieldGridr);

  std::cout << "got first grid " << dtsEfieldGridr.front() << std::endl;
  int d3 = read_profile1d(
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridZString"),
      dtsEfieldGridz);

  std::cout << "got second grid" << dtsEfieldGridz.front() << std::endl;

  int d4 = read_profile2d(
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.sheathDTS"), dtsE);
#elif EFIELD_INTERP == 2
  int nR_dtsEfield, nZ_dtsEfield;

  int d1 = read_profileNs(
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNrString"),
      cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNzString"),
      nR_dtsEfield, nZ_dtsEfield);

  sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
  sim::Array<float> dtsE(nR_dtsEfield * nZ_dtsEfield);

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
#endif
#else
  int nR_dtsEfield = 1;
  int nZ_dtsEfield = 1;
  sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
  sim::Array<float> dtsE(nR_dtsEfield * nZ_dtsEfield);
#endif

  std::string outnameEfieldR = "EfieldR.m";
  std::string outnameEfieldZ = "EfieldZ.m";
  std::string outnameEfieldT = "EfieldT.m";
  std::string outnameMinDist = "DistToSurface.m";
  // OUTPUT2d(profiles_folder,outnameEfieldR, nR_Bfield, nZ_Bfield,
  // &Efieldr.front()); OUTPUT2d(profiles_folder,outnameEfieldZ, nR_Bfield,
  // nZ_Bfield, &Efieldz.front()); OUTPUT2d(profiles_folder,outnameEfieldT,
  // nR_Bfield, nZ_Bfield, &Efieldt.front());
  // OUTPUT2d(profiles_folder,outnameMinDist, nR_Bfield, nZ_Bfield,
  // &minDist.front());

#if SPECTROSCOPY > 0
  float netX0 = 0.0, netX1 = 0.0, netY0 = 0.0, netY1 = 0.0, netZ0 = 0.0,
        netZ1 = 0.0;
  int net_nX = 0, net_nY = 0, net_nZ = 0;
  int nBins = 0;
  int nSpec = 0;
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
#if SPECTROSCOPY < 3
  nSpec = (nBins + 1) * net_nX * net_nZ;
#else
  nSpec = (nBins + 1) * net_nX * net_nY * net_nZ;
#endif
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
#if SPECTROSCOPY < 3
  sim::Array<double> net_Bins((nBins + 1) * net_nX * net_nZ, 0.0);
  sim::Array<double> net_BinsTotal((nBins + 1) * net_nX * net_nZ, 0.0);
#else
  sim::Array<double> net_Bins((nBins + 1) * net_nX * net_nY * net_nZ);
  sim::Array<double> net_BinsTotal((nBins + 1) * net_nX * net_nY * net_nZ);
#endif

  /*
  for (int i=0; i<nBins*net_nX*net_nZ; i++)
      {
          std::cout << "i " << i << std::endl;
        net_Bins[i] = 0;
          std::cout << "net bins " << net_Bins[i] << std::endl;

      }
  */
  sim::Array<float> gridX_bins(net_nX), gridY_bins(net_nY), gridZ_bins(net_nZ);

  for (int i = 0; i < net_nX; i++) {
    gridX_bins[i] = netX0 + 1.0 / (net_nX - 1) * i * (netX1 - netX0);
  }
  for (int i = 0; i < net_nY; i++) {
    gridY_bins[i] = netY0 + 1.0 / (net_nY - 1) * i * (netY1 - netY0);
  }

  for (int i = 0; i < net_nZ; i++) {
    gridZ_bins[i] = netZ0 + i * 1.0 / (net_nZ - 1) * (netZ1 - netZ0);
  }
#endif

  // Perp DiffusionCoeff initialization - only used when Diffusion interpolator
  // is = 0
  float perpDiffusionCoeff = 0.0;
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
  // Surface model import
  int nE_sputtRefCoeff = 1, nA_sputtRefCoeff = 1;
  int nE_sputtRefDistIn = 1, nA_sputtRefDistIn = 1;
  int nE_sputtRefDistOut = 1, nA_sputtRefDistOut = 1;
  int nE_sputtRefDistOutRef = 1, nDistE_surfaceModelRef = 1;
  int nDistE_surfaceModel = 1, nDistA_surfaceModel = 1;
#if USESURFACEMODEL > 0
  std::string surfaceModelCfg = "surfaceModel.";
  std::string surfaceModelFile;
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVariable(cfg, surfaceModelCfg + "fileString", surfaceModelFile);
    nE_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
                                      surfaceModelCfg, "nEsputtRefCoeffString");
    nA_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
                                      surfaceModelCfg, "nAsputtRefCoeffString");
    nE_sputtRefDistIn =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistInString");
    nA_sputtRefDistIn =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nAsputtRefDistInString");
    nE_sputtRefDistOut =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistOutString");
    nE_sputtRefDistOutRef =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistOutStringRef");
    nA_sputtRefDistOut =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nAsputtRefDistOutString");
    nDistE_surfaceModel =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOut;
    nDistE_surfaceModelRef =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOutRef;
    nDistA_surfaceModel =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nA_sputtRefDistOut;
    std::cout << " got dimensions of surface model " << std::endl;
#if USE_MPI > 0
  }
  MPI_Bcast(&nE_sputtRefCoeff, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nA_sputtRefCoeff, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nE_sputtRefDistIn, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nA_sputtRefDistIn, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nE_sputtRefDistOut, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nE_sputtRefDistOutRef, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nA_sputtRefDistOut, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nDistE_surfaceModel, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nDistE_surfaceModelRef, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nDistA_surfaceModel, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  sim::Array<float> E_sputtRefCoeff(nE_sputtRefCoeff),
      A_sputtRefCoeff(nA_sputtRefCoeff), Elog_sputtRefCoeff(nE_sputtRefCoeff),
      energyDistGrid01(nE_sputtRefDistOut),
      energyDistGrid01Ref(nE_sputtRefDistOutRef),
      angleDistGrid01(nA_sputtRefDistOut),
      spyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
      rfyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
      E_sputtRefDistIn(nE_sputtRefDistIn), A_sputtRefDistIn(nA_sputtRefDistIn),
      Elog_sputtRefDistIn(nE_sputtRefDistIn),
      E_sputtRefDistOut(nE_sputtRefDistOut),
      E_sputtRefDistOutRef(nE_sputtRefDistOutRef),
      Aphi_sputtRefDistOut(nA_sputtRefDistOut),
      Atheta_sputtRefDistOut(nA_sputtRefDistOut),
      AphiDist_Y(nDistA_surfaceModel), AthetaDist_Y(nDistA_surfaceModel),
      EDist_Y(nDistE_surfaceModel), AphiDist_R(nDistA_surfaceModel),
      AthetaDist_R(nDistA_surfaceModel), EDist_R(nDistE_surfaceModelRef),
      AphiDist_CDF_Y(nDistA_surfaceModel),
      AthetaDist_CDF_Y(nDistA_surfaceModel), EDist_CDF_Y(nDistE_surfaceModel),
      AphiDist_CDF_R(nDistA_surfaceModel),
      AthetaDist_CDF_R(nDistA_surfaceModel),
      EDist_CDF_R(nDistE_surfaceModelRef),
      AphiDist_CDF_Y_regrid(nDistA_surfaceModel),
      AthetaDist_CDF_Y_regrid(nDistA_surfaceModel),
      EDist_CDF_Y_regrid(nDistE_surfaceModel),
      AphiDist_CDF_R_regrid(nDistA_surfaceModel),
      AthetaDist_CDF_R_regrid(nDistA_surfaceModel),
      EDist_CDF_R_regrid(nDistE_surfaceModelRef);
#if USESURFACEMODEL > 0
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefCoeff", E_sputtRefCoeff[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "A_sputtRefCoeff", A_sputtRefCoeff[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistIn", E_sputtRefDistIn[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "A_sputtRefDistIn", A_sputtRefDistIn[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistOut", E_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistOutRef", E_sputtRefDistOutRef[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "Aphi_sputtRefDistOut", Aphi_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "Atheta_sputtRefDistOut", Atheta_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "sputtYldString", spyl_surfaceModel[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "reflYldString", rfyl_surfaceModel[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "EDist_Y", EDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AphiDist_Y", AphiDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AthetaDist_Y", AthetaDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "EDist_R", EDist_R[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AphiDist_R", AphiDist_R[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AthetaDist_R", AthetaDist_R[0]);
    // for(int i=0;i<nDistE_surfaceModel;i++)
    //{
    //    std::cout << " Edist diff Y " << EDist_Y[i] << " " << EDist_R[i] <<
    //    std::endl;
    //}
    for (int i = 0; i < nE_sputtRefCoeff; i++) {
      Elog_sputtRefCoeff[i] = log10(E_sputtRefCoeff[i]);
      std::cout << " EsputtRefCoeff and Elog " << E_sputtRefCoeff[i] << " "
                << Elog_sputtRefCoeff[i] << std::endl;
    }
    for (int i = 0; i < nE_sputtRefDistIn; i++) {
      Elog_sputtRefDistIn[i] = std::log10(E_sputtRefDistIn[i]);
    }
    for (int i = 0; i < nE_sputtRefDistOut; i++) {
      energyDistGrid01[i] = i * 1.0 / nE_sputtRefDistOut;
    }
    for (int i = 0; i < nE_sputtRefDistOutRef; i++) {
      energyDistGrid01Ref[i] = i * 1.0 / nE_sputtRefDistOutRef;
    }
    for (int i = 0; i < nA_sputtRefDistOut; i++) {
      angleDistGrid01[i] = i * 1.0 / nA_sputtRefDistOut;
      // std::cout << " angleDistGrid01[i] " << angleDistGrid01[i] << std::endl;
    }
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
              EDist_Y.data(), EDist_CDF_Y.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AphiDist_Y.data(), AphiDist_CDF_Y.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_Y.data(), AthetaDist_CDF_Y.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
              EDist_R.data(), EDist_CDF_R.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AphiDist_R.data(), AphiDist_CDF_R.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_R.data(), AthetaDist_CDF_R.data());
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_R.data(), AthetaDist_CDF_R.data());
    // for(int k=0;k<nE_sputtRefDistOut;k++)
    //{
    //      std::cout << "Edist_CDF_Y " <<
    //      EDist_CDF_Y[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    //      0*nE_sputtRefDistOut+k] << std::endl;
    ////      std::cout << "cosDist_CDFR " <<
    ///EDist_CDF_R[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    ///0*nE_sputtRefDistOut+k] << std::endl;
    //}
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AphiDist_CDF_Y.data(), AphiDist_CDF_Y_regrid.data());
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AthetaDist_CDF_Y.data(), AthetaDist_CDF_Y_regrid.data());
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
                energyDistGrid01.data(), nE_sputtRefDistOut,
                E_sputtRefDistOut[nE_sputtRefDistOut - 1], EDist_CDF_Y.data(),
                EDist_CDF_Y_regrid.data());
    // std::cout << "max value " << E_sputtRefDistOut[nE_sputtRefDistOut-1] <<
    // std::endl; for(int k=0;k<60;k++)
    // {
    //     std::cout << "Edis amd cdf " << k << " " << EDist_CDF_Y[k] << " "
    //     <<EDist_CDF_Y_regrid[k] << std::endl;
    // }
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AphiDist_CDF_R.data(), AphiDist_CDF_R_regrid.data());
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AthetaDist_CDF_R.data(), AthetaDist_CDF_R_regrid.data());
    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
                energyDistGrid01Ref.data(), nE_sputtRefDistOutRef,
                E_sputtRefDistOutRef[nE_sputtRefDistOutRef - 1],
                EDist_CDF_R.data(), EDist_CDF_R_regrid.data());
    // regrid2dCDF(nE_surfaceModel,nA_surfaceModel,nEdistBins_surfaceModel,energyDistGrid01.data(),nEdistBins_surfaceModel,100.0,energyDist_CDF.data(),energyDist_CDFregrid.data());
    // for(int k=0;k<nE_sputtRefDistOut;k++)
    // {
    //       std::cout << "EDist_CDFregridY " <<
    //       EDist_CDF_Y_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    //       0*nE_sputtRefDistOut+k] << std::endl;
    ////       std::cout << "cosDist_CDFregridR " <<
    ///EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    ///0*nE_sputtRefDistOut+k] << std::endl;
    // }
    // for(int k=0;k<nA_sputtRefDistOut;k++)
    // {
    //       std::cout << "ADist_CDFregridY " << k << " " << AphiDist_Y[k]<< " "
    //       << AphiDist_CDF_Y[0*nA_sputtRefDistIn*nA_sputtRefDistOut +
    //       0*nA_sputtRefDistOut+k]<< " " <<
    //       AphiDist_CDF_Y_regrid[0*nA_sputtRefDistIn*nA_sputtRefDistOut +
    //       0*nA_sputtRefDistOut+k] << std::endl;
    ////       std::cout << "cosDist_CDFregridR " <<
    ///EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    ///0*nE_sputtRefDistOut+k] << std::endl;
    // }
    // for(int k=0;k<nA_sputtRefDistOut;k++)
    // {
    //       std::cout << "ADist_CDFregridR " << k << " " << AthetaDist_R[k]<< "
    //       " << AthetaDist_CDF_R[0*nA_sputtRefDistIn*nA_sputtRefDistOut +
    //       0*nA_sputtRefDistOut+k]<< " " <<
    //       AthetaDist_CDF_R_regrid[0*nA_sputtRefDistIn*nA_sputtRefDistOut +
    //       0*nA_sputtRefDistOut+k] << std::endl;
    ////       std::cout << "cosDist_CDFregridR " <<
    ///EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut +
    ///0*nE_sputtRefDistOut+k] << std::endl;
    // }
    // float spylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff,
    // nE_sputtRefCoeff,A_sputtRefCoeff.data(),
    //                          Elog_sputtRefCoeff.data(),spyl_surfaceModel.data());
    // float rfylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff,
    // nE_sputtRefCoeff,A_sputtRefCoeff.data(),
    //                        Elog_sputtRefCoeff.data(),rfyl_surfaceModel.data());
    float spylAInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AphiDist_CDF_Y_regrid.data());
    float spylAthetaInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AthetaDist_CDF_Y_regrid.data());
    float sputEInterpVal = interp3d(
        0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), EDist_CDF_Y_regrid.data());
    float rfylAInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AphiDist_CDF_R_regrid.data());
    float rfylAthetaInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AthetaDist_CDF_R_regrid.data());
    float rflEInterpVal = interp3d(
        0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), EDist_CDF_R_regrid.data());
    // float rflAInterpVal = interp3d (
    // 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
    //       angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data()
    //       ,ADist_CDF_R_regrid.data() );
    // float rflEInterpVal = interp3d (
    // 0.44,5.0,log10(250.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
    //         energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data()
    //         ,EDist_CDF_R_regrid.data() );
    // std::cout << "Finished surface model import " <<spylInterpVal << " " <<
    // spylAInterpVal << " " << sputEInterpVal << " "<< rfylInterpVal<< " " <<
    // rflAInterpVal << " " << rflEInterpVal <<  std::endl;
    std::cout << "Finished surface model import sputtering" << spylAInterpVal
              << " " << spylAthetaInterpVal << " " << sputEInterpVal
              << std::endl;
    std::cout << "Finished surface model import reflection" << rfylAInterpVal
              << " " << rfylAthetaInterpVal << " " << rflEInterpVal
              << std::endl;
#if USE_MPI > 0
  }
  MPI_Bcast(E_sputtRefCoeff.data(), nE_sputtRefCoeff, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(A_sputtRefCoeff.data(), nA_sputtRefCoeff, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(Elog_sputtRefCoeff.data(), nE_sputtRefCoeff, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(energyDistGrid01.data(), nE_sputtRefDistOut, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(energyDistGrid01Ref.data(), nE_sputtRefDistOutRef, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(angleDistGrid01.data(), nA_sputtRefDistOut, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(spyl_surfaceModel.data(), nE_sputtRefCoeff * nA_sputtRefCoeff,
            MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(rfyl_surfaceModel.data(), nE_sputtRefCoeff * nA_sputtRefCoeff,
            MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(E_sputtRefDistIn.data(), nE_sputtRefDistIn, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(A_sputtRefDistIn.data(), nA_sputtRefDistIn, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(Elog_sputtRefDistIn.data(), nE_sputtRefDistIn, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(E_sputtRefDistOut.data(), nE_sputtRefDistOut, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(E_sputtRefDistOutRef.data(), nE_sputtRefDistOutRef, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(Aphi_sputtRefDistOut.data(), nA_sputtRefDistOut, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(Atheta_sputtRefDistOut.data(), nA_sputtRefDistOut, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_Y.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_Y.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_Y.data(), nDistE_surfaceModel, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_R.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_R.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_R.data(), nDistE_surfaceModelRef, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_CDF_Y.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_CDF_Y.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_CDF_Y.data(), nDistE_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_CDF_R.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_CDF_R.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_CDF_R.data(), nDistE_surfaceModelRef, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_CDF_Y_regrid.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_CDF_Y_regrid.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_CDF_Y_regrid.data(), nDistE_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AphiDist_CDF_R_regrid.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(AthetaDist_CDF_R_regrid.data(), nDistA_surfaceModel, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(EDist_CDF_R_regrid.data(), nDistE_surfaceModelRef, MPI_FLOAT, 0,
            MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  // Particle time stepping control
  // int ionization_nDtPerApply  =
  // cfg.lookup("timeStep.ionization_nDtPerApply"); int collision_nDtPerApply  =
  // cfg.lookup("timeStep.collision_nDtPerApply");

#ifdef __CUDACC__
  cout << "Using THRUST" << endl;
#else
  cout << "Not using THRUST" << endl;
  // int nthreads, tid;
  //#pragma omp parallel private(nthreads, tid)
  //{
  //    nthreads = omp_get_num_threads();
  //      tid = omp_get_thread_num();
  //      if(tid == 0)
  //      {
  //          std::cout << "N Threads " << nthreads << std::endl;
  //      }
  //      std::cout << "Hello world" << tid << std::endl;
  //}
  //    //nthreads = omp_get_num_threads();
  //    //nthreads = 24;
  //    //std::cout << "N threads " << nthreads << std::endl;
  // thrust::counting_iterator<std::size_t> ex0(0);
  // thrust::counting_iterator<std::size_t> ex1(nthreads-1);
  //              thrust::for_each(thrust::device, ex0,ex1,
  //                               ompPrint());
#endif

  float dt;
  int nP = 0;          // cfg.lookup("impurityParticleSource.nP");
  long nParticles = 0; // nP;
  int nT = 0;

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
  std::cout << "World rank " << world_rank << " has " << nPPerRank[world_rank]
            << " starting at " << pStartIndx[world_rank] << std::endl;
  auto particleArray = new Particles(nParticles);
  // auto particleArray2 = new Particles(nParticles);

  float x, y, z, E, vtotal, vx, vy, vz, Ex, Ey, Ez, amu, Z, charge, phi, theta,
      Ex_prime, Ez_prime, theta_transform;
  if (world_rank == 0) {
    if (cfg.lookupValue("impurityParticleSource.initialConditions.impurity_amu",
                        amu) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.impurity_Z",
                        Z) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.charge",
                        charge)) {
      std::cout << "Impurity amu Z charge: " << amu << " " << Z << " " << charge
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get point source impurity initial conditions"
          << std::endl;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&amu, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Z, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&charge, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int nSourceSurfaces = 0;
#if PARTICLE_SOURCE_SPACE == 0 // Point Source
  if (world_rank == 0) {
    if (cfg.lookupValue("impurityParticleSource.initialConditions.x_start",
                        x) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.y_start",
                        y) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.z_start",
                        z)) {
      std::cout << "Impurity point source: " << x << " " << y << " " << z
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get point source impurity initial conditions"
          << std::endl;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&y, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&z, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#elif PARTICLE_SOURCE_SPACE > 0 // Material Surfaces - flux weighted source
  Config cfg_particles;
  std::string particleSourceFile;
  getVariable(cfg, "particleSource.fileString", particleSourceFile);
  std::cout << "Open particle source file " << input_path + particleSourceFile
            << std::endl;
  importLibConfig(cfg_particles, input_path + particleSourceFile);
  std::cout << "Successfully staged input and particle source file "
            << std::endl;

  Setting &particleSourceSetting = cfg_particles.lookup("particleSource");
  std::cout << "Successfully set particleSource setting " << std::endl;
  int nSourceBoundaries = 0, nSourceElements = 0;
  float sourceMaterialZ = 0.0, accumulatedLengthArea = 0.0,
        sourceSampleResolution = 0.0;
  if (cfg_particles.lookupValue("particleSource.materialZ", sourceMaterialZ)) {
    std::cout << "Particle Source Material Z: " << sourceMaterialZ << std::endl;
  } else {
    std::cout << "ERROR: Could not get particle source material Z" << std::endl;
  }
  if (sourceMaterialZ > 0.0) { // count source boundaries
    for (int i = 0; i < nLines; i++) {
      if (boundaries[i].Z == sourceMaterialZ) {
        nSourceBoundaries++;
#if USE3DTETGEOM
        accumulatedLengthArea = accumulatedLengthArea + boundaries[i].area;
#else
        accumulatedLengthArea = accumulatedLengthArea + boundaries[i].length;
#endif
      }
    }
  } else {
    if (cfg_particles.lookupValue("particleSource.nSourceBoundaries",
                                  nSourceBoundaries)) {
      std::cout << "Particle Source nSourceBoundaries: " << nSourceBoundaries
                << std::endl;
    } else {
      std::cout << "ERROR: Could not get particle source nSourceBoundaries"
                << std::endl;
    }
    for (int i = 0; i < nSourceBoundaries; i++) {
#if USE3DTETGEOM
      accumulatedLengthArea =
          accumulatedLengthArea +
          boundaries[int(particleSourceSetting["surfaceIndices"][i])].area;
#else
      accumulatedLengthArea =
          accumulatedLengthArea +
          boundaries[int(particleSourceSetting["surfaceIndices"][i])].length;
#endif
    }
  }
  if (cfg_particles.lookupValue("particleSource.sourceSampleResolution",
                                sourceSampleResolution)) {
    std::cout << "Particle Source sample resolution: " << sourceSampleResolution
              << std::endl;
  } else {
    std::cout << "ERROR: Could not get particle source sample resolution"
              << std::endl;
  }
  nSourceElements = ceil(accumulatedLengthArea / sourceSampleResolution);
  std::cout << "nSourceBoundaries accumulatedLength nSourceElements "
            << nSourceBoundaries << " " << accumulatedLengthArea << " "
            << nSourceElements << std::endl;
  sim::Array<float> particleSourceSpaceCDF(nSourceElements, 0.0),
      particleSourceX(nSourceElements, 0.0),
      particleSourceY(nSourceElements, 0.0),
      particleSourceZ(nSourceElements, 0.0),
      particleSourceSpaceGrid(nSourceElements, 0.0);
  sim::Array<int> particleSourceIndices(nSourceElements, 0),
      particleSourceBoundaryIndices(nSourceBoundaries, 0);
#if PARTICLE_SOURCE_SPACE == 1
  for (int i = 0; i < nSourceBoundaries; i++) {
    particleSourceBoundaryIndices[i] =
        particleSourceSetting["surfaceIndices"][i];
  }
  int currentSegmentIndex = 0, currentBoundaryIndex = 0;
  float currentAccumulatedLengthArea = 0.0, lengthAlongBoundary = 0.0,
        bDotSurfaceNorm = 0.0;
  float parVec[3] = {0.0};
  float perpVec[3] = {0.0};
  currentBoundaryIndex = particleSourceBoundaryIndices[currentSegmentIndex];
  currentAccumulatedLengthArea =
      currentAccumulatedLengthArea + boundaries[currentBoundaryIndex].length;
  for (int i = 0; i < nSourceElements; i++) {
    if (i * sourceSampleResolution > currentAccumulatedLengthArea) {
      currentSegmentIndex++;
      currentBoundaryIndex = particleSourceBoundaryIndices[currentSegmentIndex];
      currentAccumulatedLengthArea = currentAccumulatedLengthArea +
                                     boundaries[currentBoundaryIndex].length;
    }
    particleSourceIndices[i] = currentBoundaryIndex;
    particleSourceBoundaryIndices[currentSegmentIndex] =
        particleSourceSetting["surfaceIndices"][currentSegmentIndex];
    boundaries[currentBoundaryIndex].getSurfaceParallel(parVec);
    lengthAlongBoundary =
        i * sourceSampleResolution - (currentAccumulatedLengthArea -
                                      boundaries[currentBoundaryIndex].length);
    particleSourceX[i] =
        boundaries[currentBoundaryIndex].x1 + parVec[0] * lengthAlongBoundary;
    particleSourceZ[i] =
        boundaries[currentBoundaryIndex].z1 + parVec[2] * lengthAlongBoundary;
    float localN = interp2dCombined(particleSourceX[i], 0.0, particleSourceZ[i],
                                    nR_Dens, nZ_Dens, DensGridr.data(),
                                    DensGridz.data(), ni.data());
    float localT = interp2dCombined(particleSourceX[i], 0.0, particleSourceZ[i],
                                    nR_Temp, nZ_Temp, TempGridr.data(),
                                    TempGridz.data(), ti.data());
    float localCs = std::sqrt(2 * localT * 1.602e-19 / (1.66e-27 * background_amu));
    float localBnorm[3] = {0.0};
    interp2dVector(&localBnorm[0], particleSourceX[i], 0.0, particleSourceZ[i],
                   nR_Bfield, nZ_Bfield, bfieldGridr.data(), bfieldGridz.data(),
                   br.data(), bz.data(), by.data());
    vectorNormalize(localBnorm, localBnorm);
    boundaries[currentBoundaryIndex].getSurfaceNormal(perpVec);
    bDotSurfaceNorm = std::abs(vectorDotProduct(localBnorm, perpVec));
    float localY = interp2dCombined(
        std::log10(3.0 * localT), 0.0, std::acos(bDotSurfaceNorm) * 180 / 3.14159,
        nE_surfaceModel, nA_surfaceModel, Elog_surfaceModel.data(),
        A_surfaceModel.data(), spyl_surfaceModel.data());
    localY = interp2dCombined(
        std::acos(bDotSurfaceNorm) * 180 / 3.14159, 0.0, std::log10(3.0 * localT),
        nA_surfaceModel, nE_surfaceModel, A_surfaceModel.data(),
        Elog_surfaceModel.data(), spyl_surfaceModel.data());
    std::cout << "LocalPotential localAngle localY " << 3.0 * localT << " "
              << std::acos(bDotSurfaceNorm) * 180 / 3.1415 << " " << localY
              << std::endl;
    float localFlux = localCs * localN * bDotSurfaceNorm; // dotB*surf
    std::cout << "segment boundary pos x z n t cs flux " << i << " "
              << currentBoundaryIndex << " " << particleSourceX[i] << " "
              << particleSourceZ[i] << " " << localN << " " << localT << " "
              << localCs << " " << localFlux << std::endl;
    std::cout << "bfield perpvec bDotSurf " << localBnorm[0] << " "
              << localBnorm[1] << " " << localBnorm[2] << " " << perpVec[0]
              << " " << perpVec[1] << " " << perpVec[2] << " "
              << bDotSurfaceNorm << " " << std::acos(bDotSurfaceNorm) * 180 / 3.1415
              << " " << localY << std::endl;
    if (i == 0) {
      particleSourceSpaceCDF[i] = localFlux * localY;
    } else {
      particleSourceSpaceCDF[i] =
          particleSourceSpaceCDF[i - 1] + localFlux * localY;
    }
    std::cout << "particleSourceSpaceCDF " << i << " "
              << particleSourceSpaceCDF[i] << std::endl;
  }
  for (int i = 0; i < nSourceElements; i++) {
    particleSourceSpaceCDF[i] =
        particleSourceSpaceCDF[i] / particleSourceSpaceCDF[nSourceElements - 1];
    std::cout << "particleSourceSpaceCDF " << i << " "
              << particleSourceIndices[i] << " " << particleSourceX[i] << " "
              << particleSourceZ[i] << particleSourceSpaceCDF[i] << std::endl;
  }
  std::random_device randDevice;
  //boost::random::mt19937 s0;
  s0.seed(123456);
  //boost::random::uniform_01<> dist01;
  float rand0 = 0.0;
  int lowInd = 0;
  int currentSegment = 0;
#else
#endif
#endif
#if PARTICLE_SOURCE_ENERGY == 0
  if (world_rank == 0) {
    if (cfg.lookupValue("impurityParticleSource.initialConditions.energy_eV",
                        E)) {
      std::cout << "Impurity point source E: " << E << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get point source impurity initial conditions"
          << std::endl;
    }
  }

#if USE_MPI > 0
  MPI_Bcast(&E, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#elif PARTICLE_SOURCE_ENERGY > 0
#if PARTICLE_SOURCE_ENERGY == 1
  // Create Thompson Distribution
  float surfaceBindingEnergy =
      cfg.lookup("impurityParticleSource.source_material_SurfaceBindingEnergy");
  float surfaceAlpha =
      cfg.lookup("impurityParticleSource.source_materialAlpha");
  std::cout << "surface binding energy " << surfaceBindingEnergy << std::endl;
  int nThompDistPoints = 200;
  float max_Energy = 100.0;
  sim::Array<float> ThompsonDist(nThompDistPoints),
      CumulativeDFThompson(nThompDistPoints);
  for (int i = 0; i < nThompDistPoints; i++) {
    if (surfaceAlpha > 0.0) {
      ThompsonDist[i] =
          surfaceAlpha * (surfaceAlpha - 1.0) *
          (i * max_Energy / nThompDistPoints) *
          std::pow(surfaceBindingEnergy, surfaceAlpha - 1.0) /
          std::pow((i * max_Energy / nThompDistPoints) + surfaceBindingEnergy,
              (surfaceAlpha + 1.0));
    } else {
      ThompsonDist[i] =
          (i * max_Energy / nThompDistPoints) /
          std::pow((i * max_Energy / nThompDistPoints) + surfaceBindingEnergy, 3);
    }
    if (i == 0) {
      CumulativeDFThompson[i] = ThompsonDist[i];
    } else {
      CumulativeDFThompson[i] = CumulativeDFThompson[i - 1] + ThompsonDist[i];
    }
  }
  for (int i = 0; i < nThompDistPoints; i++) {
    CumulativeDFThompson[i] =
        CumulativeDFThompson[i] / CumulativeDFThompson[nThompDistPoints - 1];
    // std::cout << "energy and CDF" << i*max_Energy/nThompDistPoints << " " <<
    // CumulativeDFThompson[i] << std::endl;
  }
#elif PARTICLE_SOURCE_ENERGY == 2
#endif
  //boost::random::mt19937 sE;
  //boost::random::uniform_01<> dist01E;
  float randE = 0.0;
  int lowIndE = 0;
#endif
#if PARTICLE_SOURCE_ANGLE == 0
  if (world_rank == 0) {
    if (cfg.lookupValue("impurityParticleSource.initialConditions.phi", phi) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.theta",
                        theta)) {
      std::cout << "Impurity point source angles phi theta: " << phi << " "
                << theta << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get point source angular initial conditions"
          << std::endl;
    }
  }
#if USE_MPI > 0
  MPI_Bcast(&phi, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&theta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  phi = phi * 3.141592653589793 / 180.0;
  theta = theta * 3.141592653589793 / 180.0;
  vtotal = std::sqrt(2.0 * E * 1.602e-19 / amu / 1.66e-27);
  vx = vtotal * std::sin(phi) * std::cos(theta);
  vy = vtotal * std::sin(phi) * std::sin(theta);
  vz = vtotal * std::cos(phi);
  if (phi == 0.0) {
    vx = 0.0;
    vy = 0.0;
    vz = vtotal;
  }
#elif PARTICLE_SOURCE_ANGLE > 0

  std::cout << "Read particle source " << std::endl;
#if PARTICLE_SOURCE_ENERGY < 2
  Config cfg_particles;
#endif
  // cfg_particles.readFile((input_path+"particleSource.cfg").c_str());
  // Setting& particleSource = cfg_particles.lookup("particleSource");
  // int nSegmentsAngle = particleSource["nSegmentsAngle"];
  // float angleSample;
  // sim::Array<float> sourceAngleSegments(nSegmentsAngle);
  // sim::Array<float> angleCDF(nSegmentsAngle);
  // for (int i=0; i<(nSegmentsAngle); i++)
  //{
  //    sourceAngleSegments[i] = particleSource["angles"][i];
  //    angleCDF[i] = particleSource["angleCDF"][i];
  //}
  std::random_device randDevice_particleA;
  std::mt19937 sA(randDevice_particleA());
  std::uniform_real_distribution<float> dist01A(0.0, 1.0);
  float randA = 0.0;
  int lowIndA = 0;
#endif
  std::cout << "Starting psourcefile import " << std::endl;
#if PARTICLE_SOURCE_FILE > 0 // File source
  libconfig::Config cfg_particles;
  vector<float> xpfile(nP), ypfile(nP), zpfile(nP), vxpfile(nP), vypfile(nP),
      vzpfile(nP);
  std::string ncParticleSourceFile;
  int nPfile = 0;
  if (world_rank == 0) {
    getVariable(cfg, "particleSource.ncFileString", ncParticleSourceFile);
    std::cout << "About to try to open NcFile ncp0 " << std::endl;
    // Return this in event of a problem.
    static const int NC_ERR = 2;
    try {
      NcFile ncp0("input/" + ncParticleSourceFile, NcFile::read);
    } catch (NcException &e) {
      e.what();
      cout << "FAILURE*************************************" << endl;
      return NC_ERR;
    }
    std::cout << "finished NcFile ncp0 starting ncp" << std::endl;
    NcFile ncp("input/" + ncParticleSourceFile, NcFile::read);
    std::cout << "getting dim nP" << std::endl;
    NcDim ps_nP(ncp.getDim("nP"));

    nPfile = ps_nP.getSize();
    xpfile.resize(nPfile);
    ypfile.resize(nPfile);
    zpfile.resize(nPfile);
    vxpfile.resize(nPfile);
    vypfile.resize(nPfile);
    vzpfile.resize(nPfile);
    // std::cout << "nPfile "<< nPfile << std::endl;
    NcVar ncp_x(ncp.getVar("x"));
    NcVar ncp_y(ncp.getVar("y"));
    NcVar ncp_z(ncp.getVar("z"));
    NcVar ncp_vx(ncp.getVar("vx"));
    NcVar ncp_vy(ncp.getVar("vy"));
    NcVar ncp_vz(ncp.getVar("vz"));
    std::cout << "got through NcVar " << std::endl;
    ncp_x.getVar(&xpfile[0]);
    ncp_y.getVar(&ypfile[0]);
    ncp_z.getVar(&zpfile[0]);
    ncp_vx.getVar(&vxpfile[0]);
    ncp_vy.getVar(&vypfile[0]);
    ncp_vz.getVar(&vzpfile[0]);
    std::cout << "defined file vectors " << std::endl;
    ncp.close();
    std::cout << "closed ncp " << std::endl;
  }
#if USE_MPI > 0
  MPI_Bcast(&nPfile, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank > 0) {
    xpfile.resize(nPfile);
    ypfile.resize(nPfile);
    zpfile.resize(nPfile);
    vxpfile.resize(nPfile);
    vypfile.resize(nPfile);
    vzpfile.resize(nPfile);
  }
  MPI_Bcast(xpfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(ypfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(zpfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vxpfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vypfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(vzpfile.data(), nPfile, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  // for(int i=0;i<nPfile;i++)
  //{
  //    std::cout << " xyz from file " << xpfile[i] << " " << ypfile[i] << " "
  //    << zpfile[i] << std::endl; std::cout << " Exyz from file " << Expfile[i]
  //    << " " << Eypfile[i] << " " << Ezpfile[i] << std::endl;
  //}
#endif
  std::cout << "particle file import done" << std::endl;
#if USE3DTETGEOM > 0
  // MPI_Bcast(&boundaries[0].area, nLines,MPI_FLOAT,0,MPI_COMM_WORLD);
#endif
  sim::Array<float> pSurfNormX(nP), pSurfNormY(nP), pSurfNormZ(nP), px(nP),
      py(nP), pz(nP), pvx(nP), pvy(nP), pvz(nP);
  int surfIndexMod = 0;
  float eVec[3] = {0.0};
  for (int i = 0; i < nP; i++) {
  //std::cout<< "setting particle " << i << std::endl;
#if PARTICLE_SOURCE_SPACE > 0 // File source
#if USE3DTETGEOM > 0
    surfIndexMod = i % nSourceSurfaces;
    float xCentroid = (boundaries[sourceElements[surfIndexMod]].x1 +
                       boundaries[sourceElements[surfIndexMod]].x2 +
                       boundaries[sourceElements[surfIndexMod]].x3) /
                      3.0;
    float yCentroid = (boundaries[sourceElements[surfIndexMod]].y1 +
                       boundaries[sourceElements[surfIndexMod]].y2 +
                       boundaries[sourceElements[surfIndexMod]].y3) /
                      3.0;
    float zCentroid = (boundaries[sourceElements[surfIndexMod]].z1 +
                       boundaries[sourceElements[surfIndexMod]].z2 +
                       boundaries[sourceElements[surfIndexMod]].z3) /
                      3.0;
    float bufferLaunch = 1.0e-4;
    x = xCentroid -
        bufferLaunch * boundaries[sourceElements[surfIndexMod]].a /
            boundaries[sourceElements[surfIndexMod]]
                .plane_norm; // boundaries[sourceElements[surfIndexMod]].x1;
    y = yCentroid -
        bufferLaunch * boundaries[sourceElements[surfIndexMod]].b /
            boundaries[sourceElements[surfIndexMod]]
                .plane_norm; // boundaries[sourceElements[surfIndexMod]].y1;
    z = zCentroid -
        bufferLaunch * boundaries[sourceElements[surfIndexMod]].c /
            boundaries[sourceElements[surfIndexMod]]
                .plane_norm; // boundaries[sourceElements[surfIndexMod]].z1;
#else
    // x = sampled
    rand0 = dist01(s0);
    float distAlongSegs =
        interp1dUnstructured(rand0, nSourceElements, accumulatedLengthArea,
                             &particleSourceSpaceCDF[0], lowInd);
    currentSegment = particleSourceIndices[lowInd];
    std::cout << "rand of " << rand0 << " puts the particle " << distAlongSegs
              << " along the segments on the boundary element "
              << currentSegment << std::endl;
    float parVec[3] = {0.0};
    boundaries[currentSegment].getSurfaceParallel(parVec);
    x = particleSourceX[lowInd] + (rand0 - particleSourceSpaceCDF[lowInd]) /
                                      (particleSourceSpaceCDF[lowInd + 1] -
                                       particleSourceSpaceCDF[lowInd]) *
                                      sourceSampleResolution * parVec[0];
    y = 0.0;
    z = particleSourceZ[lowInd] + (rand0 - particleSourceSpaceCDF[lowInd]) /
                                      (particleSourceSpaceCDF[lowInd + 1] -
                                       particleSourceSpaceCDF[lowInd]) *
                                      sourceSampleResolution * parVec[2];
    float buffer = 1e-6; // 0.0;//2e-6;
    x = x - buffer * boundaries[currentSegment].a /
                boundaries[currentSegment]
                    .plane_norm; // boundaries[sourceElements[surfIndexMod]].x1;
    z = z - buffer * boundaries[currentSegment].c /
                boundaries[currentSegment]
                    .plane_norm; // boundaries[sourceElements[surfIndexMod]].z1;
#endif
#endif
#if PARTICLE_SOURCE_ENERGY > 0
    randE = dist01E(sE);
#if PARTICLE_SOURCE_ENERGY == 1
    E = interp1dUnstructured(randE, nThompDistPoints, max_Energy,
                             &CumulativeDFThompson.front(), lowIndE);
#elif PARTICLE_SOURCE_ENERGY == 2
    float localT = interp2dCombined(x, y, z, nR_Temp, nZ_Temp, TempGridr.data(),
                                    TempGridz.data(), ti.data());
    float localBnorm[3] = {0.0};
    interp2dVector(&localBnorm[0], x, y, z, nR_Bfield, nZ_Bfield,
                   bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(),
                   by.data());
    vectorNormalize(localBnorm, localBnorm);
    boundaries[currentSegment].getSurfaceNormal(perpVec);
    bDotSurfaceNorm = std::abs(vectorDotProduct(localBnorm, perpVec));
    float localAngle = std::acos(bDotSurfaceNorm) * 180 / 3.1415;
    float sputtE =
        interp3d(randE, localAngle, std::log10(3.0 * localT),
                 nEdistBins_surfaceModel, nA_surfaceModel, nE_surfaceModel,
                 energyDistGrid01.data(), A_surfaceModel.data(),
                 Elog_surfaceModel.data(), energyDist_CDFregrid.data());
    E = sputtE;
    std::cout << "randE of " << randE << " with localAngle " << localAngle
              << " and local potential " << 3.0 * localT
              << " puts the particle energy to " << E << std::endl;
#endif
#endif
#if PARTICLE_SOURCE_ANGLE == 1 // Analytic normal incidence
    Ex = -E * boundaries[currentSegment].a /
         boundaries[currentSegment].plane_norm;
    Ey = -E * boundaries[currentSegment].b /
         boundaries[currentSegment].plane_norm;
    Ez = -E * boundaries[currentSegment].c /
         boundaries[currentSegment].plane_norm;

#elif PARTICLE_SOURCE_ANGLE > 1
    randA = dist01A(sA);
    float sputtA =
        interp3d(randA, localAngle, std::log10(3.0 * localT),
                 nAdistBins_surfaceModel, nA_surfaceModel, nE_surfaceModel,
                 cosDistGrid01.data(), A_surfaceModel.data(),
                 Elog_surfaceModel.data(), cosDist_CDFregrid.data());
    phi = sputtA * 3.141592653589793 / 180.0;
    std::cout << "sputtA and phi " << sputtA << " " << phi << std::endl;
    randA = dist01A(sA);
    theta = 2.0 * 3.141592653589793 * randA;
    std::cout << "randA and theta " << randA << " " << theta << std::endl;
    Ex = E * std::sin(phi) * std::cos(theta);
    Ey = E * std::sin(phi) * std::sin(theta);
    Ez = E * std::cos(phi);
    std::cout << "randA of " << randA << " puts the particle angle phi to "
              << phi << std::endl;
    std::cout << "E of particle " << Ex << " " << Ey << " " << Ez << " "
              << std::endl;
    std::cout << "current segment and perpVec " << currentSegment << " "
              << perpVec[0] << " " << perpVec[1] << " " << perpVec[2]
              << std::endl;
    float Ezx = std::sqrt(Ez * Ez + Ex * Ex);
    float thetaEzx = atan2(Ez, Ex);
    std::cout << "Ezx thetaEzx " << Ezx << " " << thetaEzx << std::endl;
    // positive slope equals negative upward normal
    theta_transform =
        std::acos(perpVec[2]); //-std::copysign(1.0,boundaries[currentSegment].slope_dzdx)*
    // if(perpVec[2]==0.0)
    //{
    //    if(perpVec[0] > 0.0)
    //    {
    //      theta_transform = 0.5*3.141592653589793;
    //      std::cout << "Vertical line element perpVec " << perpVec[0] << " "
    //      << perpVec[1] << " " << perpVec[2] << " " << theta_transform <<
    //      std::endl;
    //    }
    //    else if(perpVec[0] < 0.0)
    //    {
    //      theta_transform = 1.5*3.141592653589793;
    //      std::cout << "Vertical line element perpVec " << perpVec[0] << " "
    //      << perpVec[1] << " " << perpVec[2] << " " << theta_transform <<
    //      std::endl;
    //    }
    //}
    Ex = Ezx * std::cos(thetaEzx - theta_transform);
    // Ey = E*sin(phi+theta_transform)*sin(theta);
    Ez = Ezx * std::sin(thetaEzx - theta_transform);
    std::cout << "theta transform " << theta_transform << std::endl;
    eVec[0] = Ex;
    eVec[1] = Ey;
    eVec[2] = Ez;
    float EdotP = vectorDotProduct(perpVec, eVec);
    if (EdotP < 0.0) {
      std::cout << "This dot product negative " << std::endl;
      Ex = -Ex;
      Ez = -Ez;
    }
    // Ex_prime = Ex*cos(theta_transform) - Ez*sin(theta_transform);
    // Ez_prime = Ex*sin(theta_transform) + Ez*cos(theta_transform);
    // Ex = Ex_prime;
    // Ez = Ez_prime;
    std::cout << "Transformed E " << Ex << " " << Ey << " " << Ez << " "
              << std::endl;
    // particleArray->setParticle(i,x, y, z, Ex, Ey,Ez, Z, amu, charge);
#endif

#if PARTICLE_SOURCE_FILE > 0 // File source
    x = xpfile[i];
    y = ypfile[i];
    z = zpfile[i];
    vx = vxpfile[i];
    vy = vypfile[i];
    vz = vzpfile[i];
#endif
    // std::cout << "particle xyz Exyz Z amu charge " << x << " " << y << " " <<
    // z << " "
    // << vx << " " << vy << " " << vz << " " << Z << " " << amu << " " <<
    // charge << " "  << std::endl;
    particleArray->setParticleV(i, x, y, z, vx, vy, vz, Z, amu, charge);
#if PARTICLE_SOURCE_SPACE > 0
    pSurfNormX[i] =
        -boundaries[currentSegment].a / boundaries[currentSegment].plane_norm;
    pSurfNormY[i] =
        -boundaries[currentSegment].b / boundaries[currentSegment].plane_norm;
    pSurfNormZ[i] =
        -boundaries[currentSegment].c / boundaries[currentSegment].plane_norm;
#endif
    px[i] = x;
    py[i] = y;
    pz[i] = z;
    pvx[i] = vx;
    pvy[i] = vy;
    pvz[i] = vz;
  }
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    std::cout << "writing particles out file" << std::endl;
    NcFile ncFile_particles("output/particleSource.nc", NcFile::replace);
    NcDim pNP = ncFile_particles.addDim("nP", nP);
    NcVar p_surfNormx = ncFile_particles.addVar("surfNormX", ncFloat, pNP);
    NcVar p_surfNormy = ncFile_particles.addVar("surfNormY", ncFloat, pNP);
    NcVar p_surfNormz = ncFile_particles.addVar("surfNormZ", ncFloat, pNP);
    NcVar p_vx = ncFile_particles.addVar("vx", ncFloat, pNP);
    NcVar p_vy = ncFile_particles.addVar("vy", ncFloat, pNP);
    NcVar p_vz = ncFile_particles.addVar("vz", ncFloat, pNP);
    NcVar p_x = ncFile_particles.addVar("x", ncFloat, pNP);
    NcVar p_y = ncFile_particles.addVar("y", ncFloat, pNP);
    NcVar p_z = ncFile_particles.addVar("z", ncFloat, pNP);
    p_surfNormx.putVar(&pSurfNormX[0]);
    p_surfNormy.putVar(&pSurfNormY[0]);
    p_surfNormz.putVar(&pSurfNormZ[0]);
    p_vx.putVar(&pvx[0]);
    p_vy.putVar(&pvy[0]);
    p_vz.putVar(&pvz[0]);
    p_x.putVar(&px[0]);
    p_y.putVar(&py[0]);
    p_z.putVar(&pz[0]);
    ncFile_particles.close();
    std::cout << "finished writing particles out file" << std::endl;
#if USE_MPI > 0
  }
#endif

#if GEOM_TRACE > 0
  std::uniform_real_distribution<float> dist2(0, 1);
  // std::random_device rd2;
  // std::default_random_engine generator2(rd2());
  float randDevice02 = 6.52E+5;
  std::default_random_engine generatorTrace(randDevice02);
  std::cout << "Randomizing velocities to trace geometry. " << std::endl;

  for (int i = 0; i < nParticles; i++) {
    float theta_trace = dist2(generatorTrace) * 2 * 3.1415;
    float phi_trace = dist2(generatorTrace) * 3.1415;
    float mag_trace = 2e3;
    particleArray->vx[i] = mag_trace * std::cos(theta_trace) * std::sin(phi_trace);
    particleArray->vy[i] = mag_trace * std::sin(theta_trace) * std::sin(phi_trace);
    particleArray->vz[i] = mag_trace * std::cos(phi_trace);
  }
#endif

#if PARTICLE_TRACKS > 0
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
  for (int i = 0; i < world_size; i++) {
    pDisplacement[i] = pStartIndx[i] * nHistoriesPerParticle;
    pHistPerNode[i] = nPPerRank[i] * nHistoriesPerParticle;
  }
  const int *displ = &pDisplacement[0];
  const int *phpn = &pHistPerNode[0];
  std::cout << "History array length " << nHistories << std::endl;
#if USE_CUDA > 0
  sim::Array<float> positionHistoryX(nHistories);
  sim::Array<float> positionHistoryXgather(nHistories, 0.0);
  sim::Array<float> positionHistoryY(nHistories);
  sim::Array<float> positionHistoryYgather(nHistories);
  sim::Array<float> positionHistoryZ(nHistories);
  sim::Array<float> positionHistoryZgather(nHistories);
  sim::Array<float> velocityHistory(nHistories);
  sim::Array<float> velocityHistoryX(nHistories);
  sim::Array<float> velocityHistoryY(nHistories);
  sim::Array<float> velocityHistoryZ(nHistories);
  sim::Array<float> velocityHistorygather(nHistories);
  sim::Array<float> velocityHistoryXgather(nHistories);
  sim::Array<float> velocityHistoryYgather(nHistories);
  sim::Array<float> velocityHistoryZgather(nHistories);
  sim::Array<float> chargeHistory(nHistories);
  sim::Array<float> chargeHistoryGather(nHistories);
  sim::Array<float> weightHistory(nHistories);
  sim::Array<float> weightHistoryGather(nHistories);
#else
  std::vector<float> positionHistoryX(nHistories);
  std::vector<float> positionHistoryXgather(nHistories, 0.0);
  std::vector<float> positionHistoryY(nHistories);
  std::vector<float> positionHistoryYgather(nHistories);
  std::vector<float> positionHistoryZ(nHistories);
  std::vector<float> positionHistoryZgather(nHistories);
  std::vector<float> velocityHistory(nHistories);
  std::vector<float> velocityHistoryX(nHistories);
  std::vector<float> velocityHistoryY(nHistories);
  std::vector<float> velocityHistoryZ(nHistories);
  std::vector<float> velocityHistorygather(nHistories);
  std::vector<float> velocityHistoryXgather(nHistories);
  std::vector<float> velocityHistoryYgather(nHistories);
  std::vector<float> velocityHistoryZgather(nHistories);
  std::vector<float> chargeHistory(nHistories);
  std::vector<float> chargeHistoryGather(nHistories);
  std::vector<float> weightHistory(nHistories);
  std::vector<float> weightHistoryGather(nHistories);
#endif
#endif
  float *finalPosX = new float[nP];
  float *finalPosY = new float[nP];
  float *finalPosZ = new float[nP];
  float *finalVx = new float[nP];
  float *finalVy = new float[nP];
  float *finalVz = new float[nP];
  float *transitTime = new float[nP];
  float *hitWall = new float[nP];

  std::cout << "Beginning random number seeds" << std::endl;
  std::uniform_real_distribution<float> dist(0, 1e6);

#if FIXEDSEEDS == 0
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::default_random_engine generator1(rd());
  std::default_random_engine generator2(rd());
  std::default_random_engine generator3(rd());
  std::default_random_engine generator4(rd());
  std::default_random_engine generator5(rd());
  std::default_random_engine generator6(rd());
#endif

  thrust::counting_iterator<std::size_t> particleBegin(pStartIndx[world_rank]);
  thrust::counting_iterator<std::size_t> particleEnd(
      pStartIndx[world_rank] + nActiveParticlesOnRank[world_rank] - 1);
  thrust::counting_iterator<std::size_t> particleOne(1);
  auto randInitStart_clock = gitr_time::now();

#if PARTICLESEEDS > 0
#if USE_CUDA
  sim::Array<curandState> state1(nParticles);
#else
  sim::Array<std::mt19937> state1(nParticles);
#endif
#if USEIONIZATION > 0 || USERECOMBINATION > 0 || USEPERPDIFFUSION > 0 ||       \
    USECOULOMBCOLLISIONS > 0 || USESURFACEMODEL > 0
#if USE_CUDA
  // if(world_rank == 0)
  //{
  int *dev_i;
  cudaMallocManaged(&dev_i, sizeof(int));
  dev_i[0] = 0;
  std::cout << " about to do curandInit" << std::endl;
  thrust::for_each(thrust::device, particleBegin, particleEnd,
                   // thrust::for_each(thrust::device,particleBegin+
                   // world_rank*nP/world_size,particleBegin +
                   // (world_rank+1)*nP/world_size-10,
                   // curandInitialize(&state1[0],randDeviceInit,0));
                   curandInitialize(&state1.front(), 0));
  std::cout << " finished curandInit" << std::endl;
  // curandInitialize cuIn(0);
  // cuIn(0);
  //}
  //}
#else
  std::random_device randDeviceInit;
  // thrust::for_each(thrust::device,particleBegin+
  // world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,
  //                     curandInitialize(&state1[0],randDeviceInit,0));
  // std::mt19937 s0(randDeviceInit);
  for (int i = world_rank * nP / world_size;
       i < (world_rank + 1) * nP / world_size; i++) {
    std::mt19937 s0(randDeviceInit());
    state1[i] = s0;
  }
#endif
#if USE_CUDA
  cudaDeviceSynchronize();
#endif
#endif
#endif
  auto randInitEnd_clock = gitr_time::now();
  std::chrono::duration<float> fsRandInit = randInitEnd_clock - randInitStart_clock;
  printf(
      "Random Number Initialize time for node %i          is %6.3f (secs) \n",
      world_rank, fsRandInit.count());

  float moveTime = 0.0;
  float geomCheckTime = 0.0;
  float ionizTime = 0.0;
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
      &closeGeom_sheath.front());
  geometry_check geometry_check0(
      particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
      nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
      n_closeGeomElements.data(), &closeGeomGridr.front(),
      &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
      nEdist, E0dist, Edist, nAdist, A0dist, Adist);
#if USE_SORT > 0
  sortParticles sort0(particleArray, nP, 0.001, dev_tt, 10000,
                      pStartIndx.data(), nActiveParticlesOnRank.data(),
                      world_rank, &state1.front());
#endif
#if SPECTROSCOPY > 0
  spec_bin spec_bin0(particleArray, nBins, net_nX, net_nY, net_nZ,
                     &gridX_bins.front(), &gridY_bins.front(),
                     &gridZ_bins.front(), &net_Bins.front(), dt);
#endif
#if USEIONIZATION > 0
  ionize ionize0(
      particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
      &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
      &rateCoeff_Ionization.front());
#endif
#if USERECOMBINATION > 0
  recombine recombine0(
      particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesRecombine,
      nDensitiesRecombine, gridTemperature_Recombination.data(),
      gridDensity_Recombination.data(), rateCoeff_Recombination.data());
#endif
#if USEPERPDIFFUSION > 0
  crossFieldDiffusion crossFieldDiffusion0(
      particleArray, dt, &state1.front(), perpDiffusionCoeff, nR_Bfield,
      nZ_Bfield, bfieldGridr.data(), &bfieldGridz.front(), &br.front(),
      &bz.front(), &by.front());
#endif
#if USECOULOMBCOLLISIONS > 0
  coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, bfieldGridr.data(),
      &bfieldGridz.front(), &br.front(), &bz.front(), &by.front());

#endif
#if USETHERMALFORCE > 0
  thermalForce thermalForce0(
      particleArray, dt, background_amu, nR_gradT, nZ_gradT, gradTGridr.data(),
      gradTGridz.data(), gradTiR.data(), gradTiZ.data(), gradTiY.data(),
      gradTeR.data(), gradTeZ.data(), gradTeY.data(), nR_Bfield, nZ_Bfield,
      bfieldGridr.data(), &bfieldGridz.front(), &br.front(), &bz.front(),
      &by.front());
#endif

#if USESURFACEMODEL > 0
  reflection reflection0(
      particleArray, dt, &state1.front(), nLines, &boundaries[0], surfaces,
      nE_sputtRefCoeff, nA_sputtRefCoeff, A_sputtRefCoeff.data(),
      Elog_sputtRefCoeff.data(), spyl_surfaceModel.data(),
      rfyl_surfaceModel.data(), nE_sputtRefDistOut, nE_sputtRefDistOutRef,
      nA_sputtRefDistOut, nE_sputtRefDistIn, nA_sputtRefDistIn,
      Elog_sputtRefDistIn.data(), A_sputtRefDistIn.data(),
      E_sputtRefDistOut.data(), E_sputtRefDistOutRef.data(),
      Aphi_sputtRefDistOut.data(), energyDistGrid01.data(),
      energyDistGrid01Ref.data(), angleDistGrid01.data(),
      EDist_CDF_Y_regrid.data(), AphiDist_CDF_Y_regrid.data(),
      EDist_CDF_R_regrid.data(), AphiDist_CDF_R_regrid.data(), nEdist, E0dist,
      Edist, nAdist, A0dist, Adist);
#endif

#if PARTICLE_TRACKS > 0
  history history0(particleArray, dev_tt, nT, subSampleFac, nP,
                   &positionHistoryX.front(), &positionHistoryY.front(),
                   &positionHistoryZ.front(), &velocityHistory.front(),
                   &velocityHistoryX.front(), &velocityHistoryY.front(),
                   &velocityHistoryZ.front(), &chargeHistory.front(),
                   &weightHistory.front());
#endif
#if FORCE_EVAL > 0
  if (world_rank == 0) {
    int nR_force, nZ_force;
    float forceX0, forceX1, forceZ0, forceZ1, testEnergy;
    std::string forceCfg = "forceEvaluation.";

    getVariable(cfg, forceCfg + "nR", nR_force);
    getVariable(cfg, forceCfg + "nZ", nZ_force);
    std::vector<float> forceR(nR_force, 0.0), forceZ(nZ_force, 0.0);
    std::vector<float> tIon(nR_force * nZ_force, 0.0),
        tRecomb(nR_force * nZ_force, 0.0);
    std::vector<float> dvEr(nR_force * nZ_force, 0.0),
        dvEz(nR_force * nZ_force, 0.0), dvEt(nR_force * nZ_force, 0.0);
    std::vector<float> dvBr(nR_force * nZ_force, 0.0),
        dvBz(nR_force * nZ_force, 0.0), dvBt(nR_force * nZ_force, 0.0);
    std::vector<float> dvCollr(nR_force * nZ_force, 0.0),
        dvCollz(nR_force * nZ_force, 0.0), dvCollt(nR_force * nZ_force, 0.0);
    std::vector<float> dvITGr(nR_force * nZ_force, 0.0),
        dvITGz(nR_force * nZ_force, 0.0), dvITGt(nR_force * nZ_force, 0.0);
    std::vector<float> dvETGr(nR_force * nZ_force, 0.0),
        dvETGz(nR_force * nZ_force, 0.0), dvETGt(nR_force * nZ_force, 0.0);
    getVariable(cfg, forceCfg + "X0", forceX0);
    getVariable(cfg, forceCfg + "X1", forceX1);
    getVariable(cfg, forceCfg + "Z0", forceZ0);
    getVariable(cfg, forceCfg + "Z1", forceZ1);
    getVariable(cfg, forceCfg + "particleEnergy", testEnergy);
    for (int i = 0; i < nR_force; i++) {
      forceR[i] = forceX0 + (forceX1 - forceX0) * i / (nR_force - 1);
    }
    for (int i = 0; i < nZ_force; i++) {
      forceZ[i] = forceZ0 + (forceZ1 - forceZ0) * i / (nZ_force - 1);
    }
    float Btotal = 0.0;
    for (int i = 0; i < nR_force; i++) {
      for (int j = 0; j < nZ_force; j++) {
        interp2dVector(&Btest[0], forceR[i], 0.0, forceZ[j], nR_Bfield,
                       nZ_Bfield, bfieldGridr.data(), bfieldGridz.data(),
                       br.data(), bz.data(), by.data());
        Btotal = vectorNorm(Btest);
        // std::cout << "node " << world_rank << "Bfield at  "<< forceR[i] << "
        // " << forceZ[j]<< " " << Btest[0] << " " << Btest[1] <<
        float testTi =
            interp2dCombined(0.0, 0.1, 0.0, nR_Temp, nZ_Temp, TempGridr.data(),
                             TempGridz.data(), ti.data());
        std::cout << "Finished Temperature import " << testVec << std::endl;
        //" " << Btest[2] << " " << Btotal << std::endl;
        particleArray->setParticle(0, forceR[i], 0.0, forceZ[j], testTi, 0.0,
                                   0.0, Z, amu, charge + 1.0);
        move_boris0(0);
#if USEIONIZATION > 0
        ionize0(0);
#endif
#if USERECOMBINATION > 0
        recombine0(0);
#endif
#if USECOULOMBCOLLISIONS > 0
        coulombCollisions0(0);
#endif
#if USETHERMALFORCE > 0
        thermalForce0(0);
#endif
        dvEr[j * nR_force + i] = move_boris0.electricForce[0];
        dvEz[j * nR_force + i] = move_boris0.electricForce[2];
        dvEt[j * nR_force + i] = move_boris0.electricForce[1];
        dvBr[j * nR_force + i] = move_boris0.magneticForce[0];
        dvBz[j * nR_force + i] = move_boris0.magneticForce[2];
        dvBt[j * nR_force + i] = move_boris0.magneticForce[1];
#if USEIONIZATION > 0
        tIon[j * nR_force + i] = ionize0.tion;
#endif
#if USERECOMBINATION > 0
        tRecomb[j * nR_force + i] = recombine0.tion;
#endif
#if USECOULOMBCOLLISIONS > 0
        dvCollr[j * nR_force + i] = coulombCollisions0.dv[0];
        dvCollz[j * nR_force + i] = coulombCollisions0.dv[2];
        dvCollt[j * nR_force + i] = coulombCollisions0.dv[1];
#endif
#if USETHERMALFORCE > 0
        dvITGr[j * nR_force + i] = thermalForce0.dv_ITGx;
        dvITGz[j * nR_force + i] = thermalForce0.dv_ITGz;
        dvITGt[j * nR_force + i] = thermalForce0.dv_ITGy;
        dvETGr[j * nR_force + i] = thermalForce0.dv_ETGx;
        dvETGz[j * nR_force + i] = thermalForce0.dv_ETGz;
        dvETGt[j * nR_force + i] = thermalForce0.dv_ETGy;
#endif
      }
    }
    std::cout << " about to write ncFile_forces " << std::endl;
    NcFile ncFile_force("output/forces.nc", NcFile::replace);
    NcDim nc_nRf = ncFile_force.addDim("nR", nR_force);
    NcDim nc_nZf = ncFile_force.addDim("nZ", nZ_force);
    vector<NcDim> forceDims;
    forceDims.push_back(nc_nZf);
    forceDims.push_back(nc_nRf);
    NcVar forceRf = ncFile_force.addVar("r", ncFloat, nc_nRf);
    NcVar forceZf = ncFile_force.addVar("z", ncFloat, nc_nZf);
    NcVar nction = ncFile_force.addVar("tIon", ncFloat, forceDims);
    NcVar nctrec = ncFile_force.addVar("tRec", ncFloat, forceDims);
    NcVar dvErf = ncFile_force.addVar("dvEr", ncFloat, forceDims);
    NcVar dvEzf = ncFile_force.addVar("dvEz", ncFloat, forceDims);
    NcVar dvEtf = ncFile_force.addVar("dvEt", ncFloat, forceDims);
    NcVar dvBrf = ncFile_force.addVar("dvBr", ncFloat, forceDims);
    NcVar dvBzf = ncFile_force.addVar("dvBz", ncFloat, forceDims);
    NcVar dvBtf = ncFile_force.addVar("dvBt", ncFloat, forceDims);
    NcVar dvCollrf = ncFile_force.addVar("dvCollr", ncFloat, forceDims);
    NcVar dvCollzf = ncFile_force.addVar("dvCollz", ncFloat, forceDims);
    NcVar dvColltf = ncFile_force.addVar("dvCollt", ncFloat, forceDims);
    NcVar dvITGrf = ncFile_force.addVar("dvITGr", ncFloat, forceDims);
    NcVar dvITGzf = ncFile_force.addVar("dvITGz", ncFloat, forceDims);
    NcVar dvITGtf = ncFile_force.addVar("dvITGt", ncFloat, forceDims);
    NcVar dvETGrf = ncFile_force.addVar("dvETGr", ncFloat, forceDims);
    NcVar dvETGzf = ncFile_force.addVar("dvETGz", ncFloat, forceDims);
    NcVar dvETGtf = ncFile_force.addVar("dvETGt", ncFloat, forceDims);
    forceRf.putVar(&forceR[0]);
    forceZf.putVar(&forceZ[0]);
    nction.putVar(&tIon[0]);
    nctrec.putVar(&tRecomb[0]);
    dvErf.putVar(&dvEr[0]);
    dvEzf.putVar(&dvEz[0]);
    dvEtf.putVar(&dvEt[0]);
    dvBrf.putVar(&dvBr[0]);
    dvBzf.putVar(&dvBz[0]);
    dvBtf.putVar(&dvBt[0]);
    dvCollrf.putVar(&dvCollr[0]);
    dvCollzf.putVar(&dvCollz[0]);
    dvColltf.putVar(&dvCollt[0]);
    dvITGrf.putVar(&dvITGr[0]);
    dvITGzf.putVar(&dvITGz[0]);
    dvITGtf.putVar(&dvITGt[0]);
    dvETGrf.putVar(&dvETGr[0]);
    dvETGzf.putVar(&dvETGz[0]);
    dvETGtf.putVar(&dvETGt[0]);
    ncFile_force.close();
    particleArray->setParticleV(0, px[0], py[0], pz[0], pvx[0], pvy[0], pvz[0],
                                Z, amu, charge);
  }
#endif

  auto start_clock = gitr_time::now();
  std::chrono::duration<float> fs1 = start_clock - gitr_start_clock;
  printf("Initialize time for node %i          is %6.3f (secs) \n", world_rank,
         fs1.count());
  float testFlowVec[3] = {0.0f};
#if USEFIELDALIGNEDVALUES > 0
  interpFieldAlignedVector(&testFlowVec[0], 1.4981, 0.0, -1.2408, nR_flowV,
                           nZ_flowV, flowVGridr.data(), flowVGridz.data(),
                           flowVr.data(), flowVz.data(), flowVt.data(),
                           nR_Bfield, nZ_Bfield, bfieldGridr.data(),
                           bfieldGridz.data(), br.data(), bz.data(), by.data());
#else
  interp2dVector(&testFlowVec[0], 1.4981, 0.0, -1.2408, nR_flowV, nZ_flowV,
                 flowVGridr.data(), flowVGridz.data(), flowVr.data(),
                 flowVz.data(), flowVt.data());
#endif

  float leakZ = 0.0;
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
  // int nDevices=0;
  // nDevices = omp_get_num_threads();
  //    unsigned int cpu_thread_id = omp_get_thread_num();
  //    unsigned int num_cpu_threads = omp_get_num_threads();
  //    printf("Number of CPU threads %d (ID %d)\n", cpu_thread_id,
  //    num_cpu_threads);
#if USE_MPI > 0
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    sim::Array<int> tmpInt(1, 1), tmpInt2(1, 1);
    // int nN=10000;
    // thrust::host_vector<int> h_vec(nN);
    // thrust::generate(h_vec.begin(), h_vec.end(), rand);
    //// transfer data to the device
    // thrust::device_vector<int> d_vec = h_vec;
    // float *d_vec2;
    // cudaMallocManaged(&d_vec2, 1000*sizeof(float));
    // std::cout << "created d_vec and cmalloc, starting init " << std::endl;
    // for(int k=0;k<1000;k++)
    //{   //std::cout << "k " << k << std::endl;
    //    d_vec2[k] = 1.0f;
    //}
    // for(int k=0;k<1000;k++)
    //{   //std::cout << "k " << k << std::endl;
    //    //d_vec2[k] = 1.0f;
    // thrust::sort(thrust::device,d_vec.begin()+world_rank*nN/world_size,
    // d_vec.begin()+ (world_rank+1)*nN/world_size-1); // sort data on the device
    //}
    //// transfer data back to host
    // thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
#ifdef __CUDACC__
    cudaDeviceSynchronize();
#endif
    for (tt; tt < nT; tt++) {
      // dev_tt[0] = tt;
#if USE_SORT > 0
      thrust::for_each(thrust::device, tmpInt.begin(), tmpInt.end(), sort0);
#ifdef __CUDACC__
      cudaDeviceSynchronize();
#endif
#endif

#if PARTICLE_TRACKS > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd, history0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif
      // std::cout << " world rank pstart nactive " << world_rank << " " <<
      // pStartIndx[world_rank] << "  " << nActiveParticlesOnRank[world_rank] <<
      // std::endl; thrust::for_each(thrust::device,particleBegin,particleOne,
      //     test_routinePp(particleArray));

      thrust::for_each(thrust::device, particleBegin, particleEnd, move_boris0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       geometry_check0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif

#if SPECTROSCOPY > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd, spec_bin0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USEIONIZATION > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd, ionize0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USERECOMBINATION > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd, recombine0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USEPERPDIFFUSION > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       crossFieldDiffusion0);
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       geometry_check0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif

#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USECOULOMBCOLLISIONS > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       coulombCollisions0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USETHERMALFORCE > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd,
                       thermalForce0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif

#if USESURFACEMODEL > 0
      thrust::for_each(thrust::device, particleBegin, particleEnd, reflection0);
#ifdef __CUDACC__
      // cudaThreadSynchronize();
#endif
#endif
    }
#if PARTICLE_TRACKS > 0
    tt = nT;
    // dev_tt[0] = tt;
    // std::cout << " tt for final history " << tt << std::endl;
    thrust::for_each(thrust::device, particleBegin, particleEnd, history0);
#endif

  // Ensure that all time step loop GPU kernels are complete before proceeding
#ifdef __CUDACC__
  cudaDeviceSynchronize();
#endif

  auto finish_clock = gitr_time::now();
  std::chrono::duration<float> fs = finish_clock - start_clock;
  printf("Time taken          is %6.3f (secs) \n", fs.count());
  printf("Time taken per step is %6.3f (secs) \n", fs.count() / (float)nT);
  // for(int i=0; i<nP;i++)
  //{
  //    std::cout << "Particle test value r1: " << i << " " <<
  //    particleArray->test[i] << std::endl;
  //}

  /*
for(int i=0; i<nP ; i++)
{
  std::cout << "particle " << i << " first rnd# " <<
      particleArray->test[i] << " and x " << particleArray->xprevious[i] <<
       " hitwall " << particleArray->hitWall[i] <<
       " trans " << particleArray->transitTime[i] << std::endl;
}
*/
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
  sim::Array<float> xGather(nP, 0.0);
  sim::Array<float> test0Gather(nP, 0.0);
  sim::Array<float> test1Gather(nP, 0.0);
  sim::Array<float> yGather(nP, 0.0);
  sim::Array<float> zGather(nP, 0.0);
  sim::Array<float> vGather(nP, 0.0);
  sim::Array<float> vxGather(nP, 0.0);
  sim::Array<float> vyGather(nP, 0.0);
  sim::Array<float> vzGather(nP, 0.0);
  sim::Array<float> hitWallGather(nP, 0.0);
  sim::Array<float> weightGather(nP, 0.0);
  sim::Array<float> chargeGather(nP, 0.0);
  sim::Array<float> firstIonizationTGather(nP, 0.0);
  sim::Array<float> firstIonizationZGather(nP, 0.0);
  sim::Array<int> hasLeakedGather(nP, 0);
  // float *x_gather = NULL;
  // if (world_rank == 0) {
  //      x_gather = malloc(sizeof(float)*nP);
  //}
  std::cout << "Reached MPI barrier for gather" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Gather(&particleArray->x[world_rank * nP / world_size], nP / world_size,
             MPI_FLOAT, &xGather[0], nP / world_size, MPI_FLOAT, 0,
             MPI_COMM_WORLD);
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
#if PARTICLE_TRACKS > 0

  std::vector<float> exampleArray(4, 0.0);
  std::vector<float> exampleArrayGather(4, 0.0);
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

  // MPI_Gatherv(&exampleArray[exDispl[world_rank]],2,MPI_FLOAT,&exampleArrayGather[0],excount,exdispl,MPI_FLOAT,0,MPI_COMM_WORLD);

  // for(int i=0;i<4;i++)
  //{
  //  std::cout << "rank " << world_rank << " val " << exampleArrayGather[i] <<
  //  std::endl;
  //}

  MPI_Barrier(MPI_COMM_WORLD);

  // for(int
  // i=pDisplacement[world_rank];i<pDisplacement[world_rank]+pHistPerNode[world_rank];i++)
  //{
  //  std::cout << "Rank i "<< i << " "  << world_rank << "z " <<
  //  positionHistoryZ[i] << std::endl;
  //}
  // std::cout << "starting particle tracks gather "<< world_rank<< " pstart "<<
  // pStartIndx[world_rank] << "nhist " << nHistoriesPerParticle << std::endl;
  // std::cout << "start gather 2 "<< world_rank<< " nppr "<<
  // nPPerRank[world_rank] << "nhist " << nHistoriesPerParticle << std::endl;
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
// if(world_rank ==0)
//{
// for(int i=0;i<401;i++)
//{
//  std::cout << "Rank " << world_rank << "z " << positionHistoryZgather[i] <<
//  nSpec,
#endif

#if SPECTROSCOPY > 0
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&net_Bins[0], &net_BinsTotal[0], nSpec, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#if (USESURFACEMODEL > 0 || FLUX_EA > 0)
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
#endif
#endif
  if (world_rank == 0) {
    auto MPIfinish_clock = gitr_time::now();
    std::chrono::duration<float> fsmpi = MPIfinish_clock - finish_clock;
    printf("Time taken for mpi reduction          is %6.3f (secs) \n",
           fsmpi.count());
  }
  //    tmp202 =  particleArray->vx[0];
  // std::cout << "memory access hitwall "
  //<< particleArray->xprevious[0] << std::endl;
  // std::cout << "transit time counting " << std::endl;
#if USE_MPI > 0
  if (world_rank == 0) {
#endif
    int totalHitWall = 0;
    for (int i = 0; i < nP; i++) {
      if (particleArray->hitWall[i] > 0.0)
        totalHitWall++;
    }
#if USE3DTETGEOM > 0
    float meanTransitTime0 = 0.0;
    /*
    for (int i=0; i<nP; i++)
    {
        std::cout << "loop " << i << std::endl;
        if(particleArray->hitWall[i] == 1.0)
        {
            meanTransitTime0 = meanTransitTime0 + particleArray->transitTime[i];
        }
    }
    */
    meanTransitTime0 = meanTransitTime0 / nP;
    int max_boundary = 0;
    float max_impacts = 0.0;
    int max_boundary1 = 0;
    float max_impacts1 = 0.0;
    float *impacts = new float[nLines];
    float *xOut = new float[nP];
    float *redeposit = new float[nLines];
    float *startingParticles = new float[nLines];
    float *surfZ = new float[nLines];
    // int nA = 90;
    // int nE = 1000;
    // float* impactEnergy = new float[nLines*nA*nE];
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

/*
sim::Array<float> tally00(nLines,0);
for (int j=0; j<nP; j++)
{
    tally00[particleArray->wallHit[j]] = tally00[particleArray->wallHit[j]] + 1;
}

std::cout << "bound 164p " << tally00[164] << std::endl;
std::cout << "bound 255p " << tally00[255] << std::endl;

std::cout << "bound 164 " << boundaries[164].impacts << std::endl;
std::cout << "bound 255 " << boundaries[255].impacts << std::endl;
*/
#else
  float *impacts = new float[nLines];
  float *startingParticles = new float[nLines];
  float *surfZ = new float[nLines];
  // float* impactEnergy = new float[nLines*1000];
  for (int i = 0; i < nLines; i++) {
    impacts[i] = boundaries[i].impacts;
    startingParticles[i] = boundaries[i].startingParticles;
    surfZ[i] = boundaries[i].Z;
  }
#endif
    // add initial particle erosion to surface counting
    int closestBoundaryIndex = 0;
    int surfIndex = 0;
    float minDistance = 0.0;
    float thisE[3] = {0.0f};
    for (int j = 0; j < nP; j++) {
      minDistance =
          getE(px[j], py[j], pz[j], thisE, boundaries.data(), nLines,
               nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
               n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
               &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
               &closeGeom_sheath.front(), closestBoundaryIndex);
      // std::cout << "Starting surf minDistance and closestBoundaryIndex " <<
      // closestBoundaryIndex << " " <<
      //    minDistance << std::endl;
      // std::cout << "Particle starting z and boundary z " << pz[j] << " " <<
      //    boundaries[closestBoundaryIndex].z1 << std::endl;
      if (boundaries[closestBoundaryIndex].Z > 0.0) {
        // std::cout << "Starting surfNumber and Z " <<
        // boundaries[closestBoundaryIndex].surfaceNumber << " " <<
        //    boundaries[closestBoundaryIndex].Z << std::endl;
        surfIndex = boundaries[closestBoundaryIndex].surfaceNumber;
        grossErosion[surfIndex] = grossErosion[surfIndex] + 1.0;
      }
    }
    //#if PARTICLE_SOURCE == 1
    // int ring1 = 0;
    // int ring2 = 0;
    // int noWall = 0;
    // float meanTransitTime = 0.0;
    //
    // for(int i=0; i<nP ; i++)
    //{
    //	if(particleArray->wallIndex[i] == boundaryIndex_ImpurityLaunch[0])
    //	{
    //		ring1++;
    //	}
    //	else if(particleArray->wallIndex[i] == boundaryIndex_ImpurityLaunch[1])
    //	{
    //		ring2++;
    //	}
    //
    //	if(particleArray->wallIndex[i] == 0)
    //	{
    //		noWall++;
    //	}
    //
    //	meanTransitTime = meanTransitTime + particleArray->transitTime[i];
    //
    //}
    // meanTransitTime = meanTransitTime/(nP-noWall);
    // std::cout << "Number of impurity particles deposited on ring 1 " << ring1
    // << std::endl; std::cout << "Number of impurity particles deposited on ring
    // 2 " << ring2 << std::endl; std::cout << "Number of impurity particles not
    // deposited " << noWall << std::endl; std::cout << "Mean transit time of
    // deposited particles " << meanTransitTime << std::endl; #endif
    ofstream outfile2;
    outfile2.open("output/positions.m");
    for (int i = 1; i < nP + 1; i++) {
      outfile2 << "Pos( " << i << ",:) = [ ";
      outfile2 << particleArray->x[i - 1] << " " << particleArray->y[i - 1]
               << " " << particleArray->z[i - 1] << " ];" << std::endl;
    }
    outfile2.close();
    // Write netCDF output for positions
    NcFile ncFile0("output/positions.nc", NcFile::replace);
    NcDim nc_nP0 = ncFile0.addDim("nP", nP);
    vector<NcDim> dims0;
    dims0.push_back(nc_nP0);

    NcVar nc_x0 = ncFile0.addVar("x", ncFloat, dims0);
    NcVar nc_y0 = ncFile0.addVar("y", ncFloat, dims0);
    NcVar nc_z0 = ncFile0.addVar("z", ncFloat, dims0);
    NcVar nc_vx0 = ncFile0.addVar("vx", ncFloat, dims0);
    NcVar nc_vy0 = ncFile0.addVar("vy", ncFloat, dims0);
    NcVar nc_vz0 = ncFile0.addVar("vz", ncFloat, dims0);
    NcVar nc_trans0 = ncFile0.addVar("transitTime", ncFloat, dims0);
    NcVar nc_impact0 = ncFile0.addVar("hitWall", ncFloat, dims0);
    NcVar nc_weight0 = ncFile0.addVar("weight", ncFloat, dims0);
    NcVar nc_charge0 = ncFile0.addVar("charge", ncFloat, dims0);
    NcVar nc_leak0 = ncFile0.addVar("hasLeaked", ncInt, dims0);
#if USE_MPI > 0
    nc_x0.putVar(&xGather[0]);
    nc_y0.putVar(&yGather[0]);
    nc_z0.putVar(&zGather[0]);
    nc_vx0.putVar(&vxGather[0]);
    nc_vy0.putVar(&vyGather[0]);
    nc_vz0.putVar(&vzGather[0]);
    nc_trans0.putVar(&particleArray->transitTime[0]);
    nc_impact0.putVar(&hitWallGather[0]);
    nc_weight0.putVar(&weightGather[0]);
    nc_charge0.putVar(&chargeGather[0]);
    nc_leak0.putVar(&hasLeakedGather[0]);
#else
  nc_x0.putVar(&particleArray->xprevious[0]);
  nc_y0.putVar(&particleArray->yprevious[0]);
  nc_z0.putVar(&particleArray->zprevious[0]);
  nc_vx0.putVar(&particleArray->vx[0]);
  nc_vy0.putVar(&particleArray->vy[0]);
  nc_vz0.putVar(&particleArray->vz[0]);
  nc_trans0.putVar(&particleArray->transitTime[0]);
  nc_impact0.putVar(&particleArray->hitWall[0]);
  nc_weight0.putVar(&particleArray->weight[0]);
  nc_charge0.putVar(&particleArray->charge[0]);
  nc_leak0.putVar(&particleArray->hasLeaked[0]);
#endif
    ncFile0.close();
    // auto particleArray2 = new Particles(1);
    // std::cout << "particleArray2 z weight"<<particleArray2->z[0] << " " <<
    // particleArray2->weight[0] << std::endl;
    // particleArray2->setP(particleArray,0,0);

    // std::cout << "particleArray2 z weight"<<particleArray2->z[0] << " " <<
    // particleArray2->weight[0] << std::endl;
    // sim::Array<thrust::pair<int,float>> pair1(100);
    // sim::Array<float> weights1(100,0.0);
    // sim::Array<float> charge1(particleArray->charge);
    // charge1=particleArray->weight;
    // for(int i=0;i<nP;i++) std::cout << " charge "<< i << " "  << charge1[i]
    // << std::endl;
    ////thrust::transform(charge1.begin(),
    // for(int i=0;i<100;i++)
    //{
    //  pair1[i].first = i;
    //  pair1[i].second = 1.0*i;
    //  weights1[i] = 1.0*i;
    // std::cout << "pair "  << " " << pair1[i].first << " " << pair1[i].second
    // << std::endl;
    ////for (auto it= pair1.begin();it !=pair1.end();it++)
    ////{
    ////  //pair1[it]=.first=1;
    ////  //pair1[it].second=1.0;
    ////  std::cout << "pair " << it << " " << pair1[it].first << " " <<
    ///pair1[it].second << std::endl;
    //}
    // thrust::sort(pair1.begin(),pair1.end(),ordering());
    // for(int i=0;i<100;i++)
    //{
    // std::cout << "pair "  << " " << pair1[i].first << " " << pair1[i].second
    // << std::endl; weights1[i] = pair1[i].second;
    //}
    // sim::Array<float> weightThreshold(1,38.0);
    // sim::Array<int> lowerBoundIndex(1,0);
    // for(int i=0;i<100;i++)
    //{
    // std::cout << "weights "  << " " << weights1[i] << " " <<
    // weightThreshold[0] << std::endl;
    //}
    // thrust::lower_bound(weights1.begin(), weights1.end(),
    //                    weightThreshold.begin(),weightThreshold.end() ,
    // 		   lowerBoundIndex.begin(),thrust::less<float>());
    // std::cout << " min index " << lowerBoundIndex[0] << " " <<
    // weights1[lowerBoundIndex[0]] << std::endl; float tmpWeight = 0.0; for(int
    // i=0;i<=lowerBoundIndex[0];i++)
    //{
    // tmpWeight = weights1[i];
    // weights1[i] = pair1[100-1-i].second;
    // weights1[100-1-i] = tmpWeight;
    //}
    // for(int i=0;i<100;i++)
    //{
    // std::cout << "weights "  << " " << weights1[i] << " " <<
    // weightThreshold[0] << std::endl;
    //}
#if (USESURFACEMODEL > 0 || FLUX_EA > 0)
    std::vector<int> surfaceNumbers(nSurfaces, 0);
    int srf = 0;
    for (int i = 0; i < nLines; i++) {
      if (boundaries[i].surface) {
        surfaceNumbers[srf] = i;

        srf = srf + 1;
      }
    }
    NcFile ncFile1("output/surface.nc", NcFile::replace);
    NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
    vector<NcDim> dims1;
    dims1.push_back(nc_nLines);

    vector<NcDim> dimsSurfE;
    dimsSurfE.push_back(nc_nLines);
    NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
    NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
    dimsSurfE.push_back(nc_nAngles);
    dimsSurfE.push_back(nc_nEnergies);
    NcVar nc_grossDep = ncFile1.addVar("grossDeposition", ncFloat, nc_nLines);
    NcVar nc_grossEro = ncFile1.addVar("grossErosion", ncFloat, nc_nLines);
    NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", ncFloat, nc_nLines);
    NcVar nc_spylCounts = ncFile1.addVar("spylCounts", ncInt, nc_nLines);
    NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", ncInt, nc_nLines);
    NcVar nc_sumParticlesStrike =
        ncFile1.addVar("sumParticlesStrike", ncInt, nc_nLines);
    NcVar nc_sumWeightStrike =
        ncFile1.addVar("sumWeightStrike", ncFloat, nc_nLines);
    nc_grossDep.putVar(&grossDeposition[0]);
    nc_surfNum.putVar(&surfaceNumbers[0]);
    nc_grossEro.putVar(&grossErosion[0]);
    nc_aveSpyl.putVar(&aveSputtYld[0]);
    nc_spylCounts.putVar(&sputtYldCount[0]);
    nc_sumParticlesStrike.putVar(&sumParticlesStrike[0]);
    nc_sumWeightStrike.putVar(&sumWeightStrike[0]);
    // NcVar nc_surfImpacts = ncFile1.addVar("impacts",ncFloat,dims1);
    // NcVar nc_surfRedeposit = ncFile1.addVar("redeposit",ncFloat,dims1);
    // NcVar nc_surfStartingParticles =
    // ncFile1.addVar("startingParticles",ncFloat,dims1); NcVar nc_surfZ =
    // ncFile1.addVar("Z",ncFloat,dims1);
    NcVar nc_surfEDist = ncFile1.addVar("surfEDist", ncFloat, dimsSurfE);
    NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", ncFloat, dimsSurfE);
    NcVar nc_surfSputtDist =
        ncFile1.addVar("surfSputtDist", ncFloat, dimsSurfE);
    // nc_surfImpacts.putVar(impacts);
    //#if USE3DTETGEOM > 0
    // nc_surfRedeposit.putVar(redeposit);
    //#endif
    // nc_surfStartingParticles.putVar(startingParticles);
    // nc_surfZ.putVar(surfZ);
    nc_surfEDist.putVar(&energyDistribution[0]);
    nc_surfReflDist.putVar(&reflDistribution[0]);
    nc_surfSputtDist.putVar(&sputtDistribution[0]);
    // NcVar nc_surfEDistGrid = ncFile1.addVar("gridE",ncDouble,nc_nEnergies);
    // nc_surfEDistGrid.putVar(&surfaces->gridE[0]);
    // NcVar nc_surfADistGrid = ncFile1.addVar("gridA",ncDouble,nc_nAngles);
    // nc_surfADistGrid.putVar(&surfaces->gridA[0]);
    ncFile1.close();
#endif
#if PARTICLE_TRACKS > 0

    // Write netCDF output for histories
    netCDF::NcFile ncFile_hist("output/history.nc", NcFile::replace);
    netCDF::NcDim nc_nT = ncFile_hist.addDim("nT", nHistoriesPerParticle);
    netCDF::NcDim nc_nP = ncFile_hist.addDim("nP", nP);
    vector<NcDim> dims_hist;
    dims_hist.push_back(nc_nP);
    dims_hist.push_back(nc_nT);
    // NcDim nc_nPnT = ncFile_hist.addDim("nPnT",nP*nT/subSampleFac);
    // dims_hist.push_back(nc_nPnT);
    netCDF::NcVar nc_x = ncFile_hist.addVar("x", ncDouble, dims_hist);
    netCDF::NcVar nc_y = ncFile_hist.addVar("y", ncDouble, dims_hist);
    netCDF::NcVar nc_z = ncFile_hist.addVar("z", ncDouble, dims_hist);

    netCDF::NcVar nc_v = ncFile_hist.addVar("v", ncDouble, dims_hist);
    netCDF::NcVar nc_vx = ncFile_hist.addVar("vx", ncDouble, dims_hist);
    netCDF::NcVar nc_vy = ncFile_hist.addVar("vy", ncDouble, dims_hist);
    netCDF::NcVar nc_vz = ncFile_hist.addVar("vz", ncDouble, dims_hist);

    netCDF::NcVar nc_charge = ncFile_hist.addVar("charge", ncDouble, dims_hist);
    netCDF::NcVar nc_weight = ncFile_hist.addVar("weight", ncDouble, dims_hist);
#if USE_MPI > 0
    // if(world_rank ==0)
    //{
    // for(int i=0;i<401;i++)
    //{
    //  std::cout << "Rank " << world_rank << "z " << positionHistoryZgather[i]
    //  << std::endl;
    //}
    //}
    nc_x.putVar(&positionHistoryXgather[0]);
    nc_y.putVar(&positionHistoryYgather[0]);
    nc_z.putVar(&positionHistoryZgather[0]);

    nc_v.putVar(&velocityHistorygather[0]);
    nc_vx.putVar(&velocityHistoryXgather[0]);
    nc_vy.putVar(&velocityHistoryYgather[0]);
    nc_vz.putVar(&velocityHistoryZgather[0]);

    nc_charge.putVar(&chargeHistoryGather[0]);
    nc_weight.putVar(&weightHistoryGather[0]);
#else
    nc_x.putVar(&positionHistoryX[0]);
    nc_y.putVar(&positionHistoryY[0]);
    nc_z.putVar(&positionHistoryZ[0]);

    nc_vx.putVar(&velocityHistoryX[0]);
    nc_vy.putVar(&velocityHistoryY[0]);
    nc_vz.putVar(&velocityHistoryZ[0]);

    nc_charge.putVar(&chargeHistory[0]);
#endif
    ncFile_hist.close();
#endif
#if SPECTROSCOPY > 0
    // Write netCDF output for density data
    NcFile ncFile("output/spec.nc", NcFile::replace);
    NcDim nc_nBins = ncFile.addDim("nBins", nBins + 1);
    NcDim nc_nR = ncFile.addDim("nR", net_nX);
#if SPECTROSCOPY > 2
    NcDim nc_nY = ncFile.addDim("nY", net_nY);
#endif
    NcDim nc_nZ = ncFile.addDim("nZ", net_nZ);
    vector<NcDim> dims;
    dims.push_back(nc_nBins);
    dims.push_back(nc_nZ);
#if SPECTROSCOPY > 2
    dims.push_back(nc_nY);
#endif
    dims.push_back(nc_nR);

    NcVar nc_n = ncFile.addVar("n", ncFloat, dims);
    NcVar nc_gridR = ncFile.addVar("gridR", ncFloat, nc_nR);
    NcVar nc_gridZ = ncFile.addVar("gridZ", ncFloat, nc_nZ);
    nc_gridR.putVar(&gridX_bins[0]);
    nc_gridZ.putVar(&gridZ_bins[0]);
#if SPECTROSCOPY > 2
    NcVar nc_gridY = ncFile.addVar("gridY", ncFloat, nc_nY);
    nc_gridY.putVar(&gridY_bins[0]);
#endif
    nc_n.putVar(&net_BinsTotal[0]);
    ncFile.close();
#endif
#ifdef __CUDACC__
    cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
/*
    for(int i=0;i<100;i++)
{
    std::cout << "tests " << particleArray->test[i] << " "<<
particleArray->test0[i] <<" "<< particleArray->test1[i] << " " <<
particleArray->test2[i] << " " << particleArray->test3[i] << " " <<
particleArray->test4[i] << std::endl;
}
    for(int i=0;i<100;i++)
{
    //std::cout << "particle ionization z and t " << firstIonizationZGather[i]
<< " " << firstIonizationTGather[i] << " " << xGather[i] << " " <<
      //vxGather[i] << " " << chargeGather[i] << std::endl;
}
*/
#endif
//    for(int i=0;i<100;i++)
//{
//    std::cout << "reflected/sputtered energy " <<
//    particleArray->newVelocity[i]   << std::endl;
//}
//#endif
#if USE_MPI > 0
  }
#endif
#ifdef __CUDACC__
  cudaError_t err = cudaDeviceReset();
// cudaProfilerStop();
#endif
#if USE_MPI > 0
  // Finalize the MPI environment.
  MPI_Finalize();
#endif
  if (world_rank == 0) {
    auto gitr_finish_clock = gitr_time::now();
    std::chrono::duration<float> fstotal = gitr_finish_clock - gitr_start_clock;
    printf("Total runtime for GITR is %6.3f (secs) \n", fstotal.count());
  }
  //#endif
  return 0;
}
