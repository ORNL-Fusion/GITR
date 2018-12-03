#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "h1.cuh"
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include "boris.h"
#include "geometryCheck.h"
#include "ionize.h"
#include "recombine.h"
#include "crossFieldDiffusion.h"
#include "coulombCollisions.h"
#include "thermalForce.h"
#include "surfaceModel.h"
#include "interp2d.hpp"
#include "interpRateCoeff.hpp"
#include <algorithm>
#include <random>
#include "Particles.h"
#include "Boundary.h"
#include "BoundaryModifiable.h"
#include "curandInitialize.h"
#include "spectroscopy.h"
#include "history.h"
#include "hashGeom.h"
#include "hashGeomSheath.h"
#include "testRoutineP.h"
#include "testRoutineCuda.h"
#include "boundaryInit.h"
#include "array.h"
#include "ompPrint.h"
#if USE_BOOST
    #include <boost/timer/timer.hpp>
    #include "boost/filesystem.hpp"
#endif

#include <vector>
#include "io.hpp"
#include "testRoutine.h"
#include "testRoutineCuda.h"
#include "boundaryInit.h"
#include <netcdf>
#include "array.h"

#ifdef __CUDACC__
    #include <curand.h>
    #include <curand_kernel.h>
#endif

#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/functional.h>

using namespace std;
using namespace libconfig;

#if USE_BOOST
    using namespace boost::timer;
#endif

using namespace netCDF;
using namespace exceptions;
using namespace netCDF::exceptions;

struct test_routinePp { 
   Particles *particles;
    test_routinePp(Particles *_particles) : particles(_particles) {}

    __device__
    void operator()(std::size_t indx) const { 
    particles->x[indx] = 5.0;
}
};

int main(int argc, char **argv)
{
  typedef std::chrono::high_resolution_clock Time;
  typedef std::chrono::duration<float> fsec;
  auto GITRstart_clock = Time::now();
  int ppn=1; 
  #if USE_MPI > 0
    ppn = 4;
    //int np = 1;
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);
  #endif
  int counter;
  printf("Program Name Is: %s",argv[0]);
  if(argc==1)
    printf("\nNo Extra Command Line Argument Passed Other Than Program Name");
  if(argc>=2)
  {
    printf("\nNumber Of Arguments Passed: %d",argc);
    printf("\n----Following Are The Command Line Arguments Passed----");
    for(counter=0;counter<argc;counter++)
    {
      printf("\nargv[%d]: %s",counter,argv[counter]);
      if(std::string(argv[counter]) == "-nGPUPerNode")
      {
        if (counter + 1 < argc) 
	{ // Make sure we aren't at the end of argv!
          ppn = std::stoi(argv[counter+1]);	  
          printf("\nGITR set to use %d GPUs per node",ppn);
	} else 
	{ // Uh-oh, there was no argument to the destination option.
	   std::cerr << "--nGPUPerNode option requires one argument." << std::endl;
          return 1;
	} 
      }
    }
  }
    
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
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);
    #if USE_CUDA > 0
      cudaSetDevice(world_rank%ppn);  
    #endif
  #else
    int world_rank=0;
    int world_size=1;
  #endif
  
  //Prepare config files for import
  Config cfg,cfg_geom;
  cfg.readFile("gitrInput.cfg");
  std::cout << "staged input file " << std::endl;
  cfg_geom.readFile("gitrGeometry.cfg");
  std::cout << "staged input files " << std::endl;
  //check binary compatibility with input file
  #if CHECK_COMPATIBILITY>0
    Setting& flags = cfg.lookup("flags");
    
    // Parse and read geometry file
    std::cout << "Open geometry file" << std::endl;
    std::string geomFile; 
    getVariable(cfg,"geometry.fileString",geomFile);
    std::cout << "Open geometry file " << input_path+geomFile << std::endl; 
    importLibConfig(cfg_geom,input_path+geomFile);

    int check3 = flags["USE_BOOST"];      
    if (USE_BOOST != check3)
    { std::cout << "incompatibility in USE_BOOST between input file and binary" << std::endl;
      exit(0);
    }
    int check4 = flags["USEIONIZATION"];      
    if (USEIONIZATION != check4)
    { std::cout << "incompatibility in USEIONIZATION between input file and binary" << std::endl;
      exit(0);
    }
    
    //check binary compatibility with input file
    #if CHECK_COMPATIBILITY>0
      checkFlags(cfg); 
    #endif
  }
  // show memory usage of GPU
  #if USE_CUDA 
  if(world_rank == 0)
  {
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t    cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  
    if(cudaSuccess != cuda_status )
    {
  
       printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
       exit(1);
    }
    
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
      used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0); 
    int nDevices;
    int nThreads;
    cudaGetDeviceCount(&nDevices);
    std::cout << "number of devices gotten " << nDevices << std::endl;
    for (int i = 0; i < nDevices; i++) {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("Device Number: %d\n", i);
      printf("  Device name: %s\n", prop.name);
      printf("  Memory Clock Rate (KHz): %d\n",
                         prop.memoryClockRate);
      printf("  Memory Bus Width (bits): %d\n",
                         prop.memoryBusWidth);
      printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
                         2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
      printf("  Total number of threads: %d\n", prop.maxThreadsPerMultiProcessor);
      nThreads = prop.maxThreadsPerMultiProcessor;
      }
  }  
  #endif
  #if USE_BOOST > 0
    std::string output_folder="output";
    //Output
    boost::filesystem::path dir(output_folder); 
    if(!(boost::filesystem::exists(dir)))
    {
      std::cout<<"Doesn't Exist in main"<<std::endl;
      if (boost::filesystem::create_directory(dir))
      {
         std::cout << " Successfully Created " << std::endl;
      }
    }

    int check7 = flags["USECOULOMBCOLLISIONS"];      
    if (USECOULOMBCOLLISIONS !=check7)
    { std::cout << "incompatibility in USECOULOMBCOLLISIONS between input file and binary" << std::endl;
      exit(0);
    }
    int check8 = flags["USETHERMALFORCE"];      
    if (USETHERMALFORCE !=check8)
    {
            std::cout << "incompatibility in USETHERMALFORCE between input file and binary" << std::endl;
            exit(0);
    }
            int check9 = flags["USESURFACEMODEL"];      
    if (USESURFACEMODEL !=check9)
    {
            std::cout << "incompatibility in USESURFACEMODEL between input file and binary" << std::endl;
            exit(0);
    }
    int check10 = flags["USESHEATHEFIELD"];      
    if (USESHEATHEFIELD !=check10)
    {
            std::cout << "incompatibility in USESHEATHEFIELD between input file and binary" << std::endl;
            exit(0);
    }
    #if USE_MPI > 0
      const int nSurfaceMembers = 18;
      
      int lengthsSurface[nSurfaceMembers] = {1,1,1,1,1,1,1,1,1,nEdist*nAdist,nEdist,nAdist,nEdist*nAdist,nEdist*nAdist,nEdist*nAdist,nEdist*nAdist,nEdist*nAdist,nSurfaces*nEdist*nAdist};
      MPI_Aint offsetsSurface[nSurfaceMembers] = {offsetof(Surfaces, nSurfaces), offsetof(Surfaces, nE),offsetof(Surfaces,nA),offsetof(Surfaces,E0), offsetof(Surfaces, E),offsetof(Surfaces,A0),
      offsetof(Surfaces, A), offsetof(Surfaces, dE),offsetof(Surfaces,dA),
      offsetof(Surfaces,sumParticlesStrike), offsetof(Surfaces, gridE),offsetof(Surfaces,gridA),
      offsetof(Surfaces,sumWeightStrike), offsetof(Surfaces, grossDeposition),offsetof(Surfaces,grossErosion),
      offsetof(Surfaces,aveSputtYld), offsetof(Surfaces,sputtYldCount),offsetof(Surfaces,energyDistribution)};
      MPI_Datatype typesSurface[nSurfaceMembers] = {MPI_INT,MPI_INT,MPI_INT,MPI_FLOAT,
      MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_INT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT};
      MPI_Datatype surface_type;
      MPI_Type_create_struct(nSurfaceMembers, lengthsSurface, offsetsSurface, typesSurface,&surface_type);    
      MPI_Type_commit(&surface_type);
      MPI_Bcast(&nEdist,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nAdist,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&E0dist,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(&Edist,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(&A0dist,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(&Adist,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    #endif
  #endif
  auto surfaces = new Surfaces(nSurfaces,nEdist,nAdist);
  surfaces->setSurface(nEdist,E0dist,Edist,nAdist,A0dist,Adist);
  
  //#if USE_MPI > 0
    //Arrays used for reduction at end of sim
    sim::Array<float> grossDeposition(nSurfaces,0.0);
    sim::Array<float> grossErosion(nSurfaces,0.0);
    sim::Array<float> sumWeightStrike(nSurfaces,0.0);
    sim::Array<float> energyDistribution(nSurfaces*nEdist*nAdist,0.0);
    sim::Array<float> reflDistribution(nSurfaces*nEdist*nAdist,0.0);
    sim::Array<float> sputtDistribution(nSurfaces*nEdist*nAdist,0.0);
    sim::Array<float> aveSputtYld(nSurfaces,0.0);
    sim::Array<int> sputtYldCount(nSurfaces,0);
    sim::Array<int> sumParticlesStrike(nSurfaces,0);
  //#endif

  int nHashes = 1;
  int nR_closeGeomTotal = 0;
  int nY_closeGeomTotal = 0;
  int nZ_closeGeomTotal = 0;
  int nHashPointsTotal = 0;
  int nGeomHash = 0;
  std::string geomHashCfg = "geometry_hash.";
  std::cout << " node starting geomhash1 " << world_rank << std::endl;
  #if GEOM_HASH == 1
    if(world_rank == 0)
    {
            std::cout << "incompatibility in USEPRESHEATHFIELD between input file and binary" << std::endl;
            exit(0);
    }
    int check12 = flags["BFIELD_INTERP"];      
    if (BFIELD_INTERP !=check12)
    {
        importHashNs(cfg,input_path,nHashes,"geometry_hash",nR_closeGeom.data(), nY_closeGeom.data(),nZ_closeGeom.data(),n_closeGeomElements.data(),nR_closeGeomTotal,nY_closeGeomTotal,nZ_closeGeomTotal,nHashPoints.data(), nHashPointsTotal,nGeomHash);
      //Setting& geomHash = cfg.lookup("geometry_hash");
      //if(nHashes > 1)
      //{
      //  for(int i=0; i<nHashes;i++)
      //  {   
      //    nR_closeGeom[i] = geomHash["nR_closeGeom"][i];
      //    nZ_closeGeom[i] = geomHash["nZ_closeGeom"][i];
      //    n_closeGeomElements[i] = geomHash["n_closeGeomElements"][i];
      //    std::cout << "hash nr ny nz total " << n_closeGeomElements[i] << " " << nR_closeGeom[i]  << " " << nZ_closeGeom[i]<< std::endl;
      //  }
      //}
      //else
      //{
      //  getVariable(cfg,geomHashCfg+"nR_closeGeom",nR_closeGeom[0]);
      //  getVariable(cfg,geomHashCfg+"nZ_closeGeom",nZ_closeGeom[0]);
      //  getVariable(cfg,geomHashCfg+"n_closeGeomElements",n_closeGeomElements[0]);
      //}
      //for(int j=0;j<nHashes;j++)
      //{
      //  nGeomHash = nGeomHash + nR_closeGeom[j]*nZ_closeGeom[j]*n_closeGeomElements[j];
      //  nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[j];
      //  nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[j];
      //}
      //#if USE3DTETGEOM > 0
      //if(nHashes > 1)
      //{
      //  for(int i=0; i<nHashes;i++)
      //  {   
      //    nY_closeGeom[i] = geomHash["nY_closeGeom"][i];
      //  }
      //}
      //else
      //{
      //  getVariable(cfg,geomHashCfg+"nY_closeGeom",nY_closeGeom[0]);
      //}
      //#endif
      //nGeomHash = 0;
      //nR_closeGeomTotal = 0;
      //nY_closeGeomTotal = 0;
      //nZ_closeGeomTotal = 0;
      //nGeomHash = 0;
      //for(int j=0;j<nHashes;j++)
      //{
      //  if(nHashes > 1)
      //  {
      //    nHashPoints[j] =nR_closeGeom[j]*nY_closeGeom[j]*nZ_closeGeom[j];
      //  }
      //  else
      //  {
      //    nHashPoints[j] =nR_closeGeom[j]*nZ_closeGeom[j];
      //  } 
      //  nHashPointsTotal = nHashPointsTotal + nHashPoints[j];
      //  nGeomHash = nGeomHash + nHashPoints[j]*n_closeGeomElements[j];
      //  nR_closeGeomTotal = nR_closeGeomTotal + nR_closeGeom[j];
      //  nY_closeGeomTotal = nY_closeGeomTotal + nY_closeGeom[j];
      //  nZ_closeGeomTotal = nZ_closeGeomTotal + nZ_closeGeom[j];
      //}
      //std::cout << "hhhash nr ny nz total " << nGeomHash << " " << nR_closeGeomTotal << " " << nY_closeGeomTotal << " " << nZ_closeGeomTotal<< std::endl;
    }
    int check13 = flags["EFIELD_INTERP"];      
    if (EFIELD_INTERP !=check13)
    {
            std::cout << "incompatibility in EFIELD_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check14 = flags["PRESHEATH_INTERP"];      
    if (PRESHEATH_INTERP !=check14)
    {
            std::cout << "incompatibility in PRESHEATH_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check15 = flags["DENSITY_INTERP"];      
    if (DENSITY_INTERP !=check15)
    {
            std::cout << "incompatibility in DENSITY_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check16 = flags["TEMP_INTERP"];      
    if (TEMP_INTERP !=check16)
    {
            std::cout << "incompatibility in TEMP_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check17 = flags["FLOWV_INTERP"];      
    if (FLOWV_INTERP !=check17)
    {
            std::cout << "incompatibility in FLOWV_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check18 = flags["GRADT_INTERP"];      
    if (GRADT_INTERP !=check18)
    {
            std::cout << "incompatibility in GRADT_INTERP between input file and binary" << std::endl;
            exit(0);
    }
    int check19 = flags["ODEINT"];      
    if (ODEINT !=check19)
    {
            std::cout << "incompatibility in ODEINT between input file and binary" << std::endl;
            exit(0);
    }
    int check20 = flags["FIXEDSEEDS"];      
    if (FIXEDSEEDS !=check20)
    {
            std::cout << "incompatibility in FIXEDSEEDS between input file and binary" << std::endl;
            exit(0);
    }
    int check21 = flags["PARTICLESEEDS"];      
    if (PARTICLESEEDS !=check21)
    {
            std::cout << "incompatibility in PARTICLESEEDS between input file and binary" << std::endl;
            exit(0);
    }
    int check22 = flags["GEOM_TRACE"];      
    if (GEOM_TRACE !=check22)
    {
            std::cout << "incompatibility in GEOM_TRACE between input file and binary" << std::endl;
            exit(0);
    }
    int check23 = flags["GEOM_HASH"];      
    if (GEOM_HASH !=check23)
    {
            std::cout << "incompatibility in GEOM_HASH between input file and binary" << std::endl;
            exit(0);
    }
    int check24 = flags["GEOM_HASH_SHEATH"];      
    if (GEOM_HASH_SHEATH !=check24)
    {
            std::cout << "incompatibility in GEOM_HASH_SHEATH between input file and binary" << std::endl;
            exit(0);
    }
    int check25 = flags["PARTICLE_TRACKS"];      
    if (PARTICLE_TRACKS !=check25)
    {
            std::cout << "incompatibility in PARTICLE_TRACKS between input file and binary" << std::endl;
            exit(0);
    }
    int check26 = flags["PARTICLE_SOURCE"];      
    if (PARTICLE_SOURCE !=check26)
    {
            std::cout << "incompatibility in PARTICLE_SOURCE between input file and binary" << std::endl;
            exit(0);
    }
    int check27 = flags["SPECTROSCOPY"];      
    if (SPECTROSCOPY !=check27)
    {
            std::cout << "incompatibility in SPECTROSCOPY between input file and binary" << std::endl;
            exit(0);
    }
    int check28 = flags["USE3DTETGEOM"];      
    if (USE3DTETGEOM !=check28)
    {
            std::cout << "incompatibility in USE3DTETGEOM between input file and binary" << std::endl;
            exit(0);
    }
    int check29 = flags["USECYLSYMM"];      
    if (USECYLSYMM !=check29)
    {
            std::cout << "incompatibility in USECYLSYMM between input file and binary" << std::endl;
            exit(0);
    }
  #endif
    std::cout << "Finished checking input file" << std::endl; 
  // Background species info
  float background_Z;
  float background_amu;
  if (cfg.lookupValue("backgroundPlasmaProfiles.Z",background_Z) &&  cfg.lookupValue("backgroundPlasmaProfiles.amu",background_amu))
  {
  float background_Z = cfg.lookup("backgroundPlasmaProfiles.Z");
  float background_amu = cfg.lookup("backgroundPlasmaProfiles.amu");
  }
  else
  {
      std::cout << "config file read failed in background Plasma profiles Z or amu " << std::endl;
  }
  //Bfield initialization
  //if(cfg.lookup("backgroundPlasmaProfiles.Bfield.br"))
  //{
  #if BFIELD_INTERP == 0
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<float> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    sim::Array<float> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);
    br[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.br");
    bz[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bz");
    bt[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bt");
  #elif BFIELD_INTERP == 2
    int nR_Bfield;
    int nZ_Bfield;
    
    int b1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNzString"),nR_Bfield,nZ_Bfield);
    
    sim::Array<float> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    sim::Array<float> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);
    
    int b2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridRString"), bfieldGridr);
    
    int b3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridZString"), bfieldGridz);
    
    int b4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.radialComponentString"), br);
    
    int b5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.axialComponentString"), bz);
    
    int b6 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.toroidalComponentString"), bt);
  #endif
  
  std::string outnameBfieldR = "BfieldR.m";
  std::string outnameBfieldZ = "BfieldZ.m";
  std::string outnameBfieldT = "BfieldT.m";
  std::string outnameGridR = "gridR.m";
  std::string outnameGridZ = "gridZ.m";
  std::string profiles_folder = "profiles";
  OUTPUT1d(profiles_folder,outnameGridR, nR_Bfield, &bfieldGridr.front());
  OUTPUT1d(profiles_folder,outnameGridZ, nZ_Bfield, &bfieldGridz.front());
  OUTPUT2d(profiles_folder,outnameBfieldR, nR_Bfield, nZ_Bfield, &br.front());
  OUTPUT2d(profiles_folder,outnameBfieldZ, nR_Bfield, nZ_Bfield, &bz.front());
  OUTPUT2d(profiles_folder,outnameBfieldT, nR_Bfield, nZ_Bfield, &bt.front());
  //}
  //else
  //{
  //    std::cout << "config file read failed in background Plasma profiles Bfield " << std::endl;
  //}

  //Background Plasma Temperature Initialization    
  #if TEMP_INTERP == 0
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<float> TempGridr(nR_Temp), TempGridz(nZ_Temp);
    sim::Array<float> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);
    ti[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.ti");
    te[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.te");
  #elif TEMP_INTERP == 2
    int nR_Temp;
    int nZ_Temp;
    
    int t1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNzString"),nR_Temp,nZ_Temp);
  
    sim::Array<float> TempGridr(nR_Temp), TempGridz(nZ_Temp);
    sim::Array<float> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);
    
    int t2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridRString"), TempGridr);
    
    int t3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridZString"), TempGridz);
    
    int t4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.IonTempString"), ti);
    
    int t5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.ElectronTempString"), te);
  #endif
  std::string outnameTi = "ti.m";
  std::string outnameTe = "te.m";
  OUTPUT2d(profiles_folder,outnameTi, nR_Temp, nZ_Temp, &ti.front());
  OUTPUT2d(profiles_folder,outnameTe, nR_Temp, nZ_Temp, &te.front());

  //Background Plasma Density Initialization
  #if DENSITY_INTERP == 0
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<float> DensGridr(nR_Dens), DensGridz(nZ_Dens);
    sim::Array<float> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
    ni[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ni");
    ne[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ne");
  #elif DENSITY_INTERP == 2
    int nR_Dens;
    int nZ_Dens;
    
    int n1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Density.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Density.gridNzString"),nR_Dens,nZ_Dens);
    
    sim::Array<float> DensGridr(nR_Dens), DensGridz(nZ_Dens);
    sim::Array<float> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
    
    int n2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.gridRString"), DensGridr);
    
    int n3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.gridZString"), DensGridz);
    
    int n4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.IonDensityString"), ni);
    
    int n5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.ElectronDensityString"), ne);
  #endif
  
  std::string outnameNi = "ni.m";
  std::string outnameNe = "ne.m";
  OUTPUT2d(profiles_folder,outnameNi, nR_Dens, nZ_Dens, &ni.front());
  OUTPUT2d(profiles_folder,outnameNe, nR_Dens, nZ_Dens, &ne.front());

  //Background Plasma flow velocity initialization    
  #if FLOWV_INTERP == 0
    int nR_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<float> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
    sim::Array<float> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),
                        flowVt(nR_flowV*nZ_flowV);
    flowVr[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVr");
    flowVz[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVz");
  #elif FLOWV_INTERP == 2
    int nR_flowV;
    int nZ_flowV;
    
    int f1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNzString"),nR_flowV,
                nZ_flowV);
    
    sim::Array<float> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
    sim::Array<float> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),
                        flowVt(nR_flowV*nZ_flowV);
    
    int f2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridRString"), flowVGridr);
    
    int f3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridZString"), flowVGridz);
    
    int f4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVrString"), flowVr);
    
    int f5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVzString"), flowVz);
    
    int f6 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVtString"), flowVt);
  #endif

  std::string outnameFlowVr = "flowVr.m";
  std::string outnameFlowVz = "flowVz.m";
  std::string outnameFlowVt = "flowVt.m";
  OUTPUT2d(profiles_folder,outnameFlowVr, nR_flowV, nZ_flowV, &flowVr.front());
  OUTPUT2d(profiles_folder,outnameFlowVz, nR_flowV, nZ_flowV, &flowVz.front());
  OUTPUT2d(profiles_folder,outnameFlowVt, nR_flowV, nZ_flowV, &flowVt.front());

  //Background plasma temperature gradient field intitialization    
  #if GRADT_INTERP == 0
    int nR_gradT = 1;
    int nZ_gradT = 1;
    sim::Array<float> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
    sim::Array<float> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),
        gradTeT(nR_gradT*nZ_gradT,0.0),gradTiR(nR_gradT*nZ_gradT), 
        gradTiZ(nR_gradT*nZ_gradT),gradTiT(nR_gradT*nZ_gradT,0.0);    
    gradTeR[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeR");
    gradTeZ[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeZ");
    gradTiR[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiR");
    gradTiZ[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiZ");
  #elif GRADT_INTERP == 2
    int nR_gradT;
    int nZ_gradT;
    
    int g1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.gradT.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.gradT.gridNzString"),nR_gradT,nZ_gradT);
    
    sim::Array<float> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
    sim::Array<float> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),
        gradTeT(nR_gradT*nZ_gradT,0.0),gradTiR(nR_gradT*nZ_gradT), 
        gradTiZ(nR_gradT*nZ_gradT),gradTiT(nR_gradT*nZ_gradT,0.0);
    
    int g2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gridRString"), gradTGridr);
    
    int g3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gridZString"), gradTGridz);
    
    int g4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiRString"), gradTiR);
    
    int g5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiZString"), gradTiZ);
    
    int g6 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeRString"), gradTeR);
    
    int g7 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeZString"), gradTeZ);
  #endif

  std::string outnameGradTiR = "gradTiR.m";
  std::string outnameGradTiZ = "gradTiZ.m";
  std::string outnameGradTeR = "gradTeR.m";
  std::string outnameGradTeZ = "gradTeZ.m";
  OUTPUT2d(profiles_folder,outnameGradTiR, nR_gradT, nZ_gradT, &gradTiR.front());
  OUTPUT2d(profiles_folder,outnameGradTiZ, nR_gradT, nZ_gradT, &gradTiZ.front());
  OUTPUT2d(profiles_folder,outnameGradTeR, nR_gradT, nZ_gradT, &gradTeR.front());
  OUTPUT2d(profiles_folder,outnameGradTeZ, nR_gradT, nZ_gradT, &gradTeZ.front());


  //Initialization of ionization and recombination coefficients    
  int nCS_Ionize, nCS_Recombine;
  int i0 = read_profileNs(cfg.lookup("impurityParticleSource.ionization.fileString"),
            cfg.lookup("impurityParticleSource.ionization.nChargeStateString"),
            cfg.lookup("impurityParticleSource.recombination.nChargeStateString"),
            nCS_Ionize, nCS_Recombine);

  int nTemperaturesIonize, nDensitiesIonize;
  int i1 = read_profileNs(cfg.lookup("impurityParticleSource.ionization.fileString"),
            cfg.lookup("impurityParticleSource.ionization.DensGridString"),
            cfg.lookup("impurityParticleSource.ionization.TempGridString"),
            nDensitiesIonize,nTemperaturesIonize);

  sim::Array<float> rateCoeff_Ionization(nCS_Ionize*nTemperaturesIonize*nDensitiesIonize);
  sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
                        gridDensity_Ionization(nDensitiesIonize);

  int i2 = read_profiles(cfg.lookup("impurityParticleSource.ionization.fileString"),
        nTemperaturesIonize,nDensitiesIonize,
        cfg.lookup("impurityParticleSource.ionization.TempGridVarName"), 
        gridTemperature_Ionization,cfg.lookup("impurityParticleSource.ionization.DensGridVarName"),
        gridDensity_Ionization,
        cfg.lookup("impurityParticleSource.ionization.CoeffVarName"),
        rateCoeff_Ionization);
   
  int nTemperaturesRecombine, nDensitiesRecombine;
  int i3 = read_profileNs(cfg.lookup("impurityParticleSource.recombination.fileString"),
            cfg.lookup("impurityParticleSource.recombination.DensGridString"),
            cfg.lookup("impurityParticleSource.recombination.TempGridString"),
            nDensitiesRecombine,nTemperaturesRecombine);

  sim::Array<float> rateCoeff_Recombination(nCS_Recombine*nTemperaturesRecombine*nDensitiesRecombine);
  sim::Array<float> gridTemperature_Recombination(nTemperaturesRecombine),
                    gridDensity_Recombination(nDensitiesRecombine);

  int i4 = read_profiles(cfg.lookup("impurityParticleSource.recombination.fileString"),
             nTemperaturesRecombine,nDensitiesRecombine,
             cfg.lookup("impurityParticleSource.recombination.TempGridVarName"), 
             gridTemperature_Recombination,cfg.lookup("impurityParticleSource.recombination.DensGridVarName"),
             gridDensity_Recombination,
             cfg.lookup("impurityParticleSource.recombination.CoeffVarName"),
             rateCoeff_Recombination);


  //Geometry Definition
  Setting& geom = cfg_geom.lookup("geom");
  int nLines = geom["x1"].getLength();
  //int nMaterials = geom["nMaterials"];
  std::cout << "Number of Geometric Objects Loaded: " << nLines << std::endl;

  auto boundaryModArray = new BoundaryModifiable(nLines);
  sim::Array<Boundary> boundaries(nLines+1);

  std::string geom_outname = "geom.m";
  std::string geom_folder = "geometry";
  ofstream outfile;

  #if USE_BOOST
    //Output
    boost::filesystem::path dir(geom_folder);
    
    if(!(boost::filesystem::exists(dir)))
    {
       std::cout<<"Doesn't Exists"<<std::endl;

       if (boost::filesystem::create_directory(dir))
       {
          std::cout << " Successfully Created " << std::endl;
       }
    }
  #endif

  std::string full_path = geom_folder + "/" + geom_outname;
  outfile.open (full_path );
  #if USE3DTETGEOM > 0
  for(int i=0 ; i<nLines ; i++)
  {
     boundaries[i].x1 = geom["x1"][i];
     boundaries[i].y1 = geom["y1"][i];
     boundaries[i].z1 = geom["z1"][i];
     boundaries[i].x2 = geom["x2"][i];
     boundaries[i].y2 = geom["y2"][i];
     boundaries[i].z2 = geom["z2"][i];
     boundaries[i].x3 = geom["x3"][i];
     boundaries[i].y3 = geom["y3"][i];
     boundaries[i].z3 = geom["z3"][i];
     boundaries[i].Z = geom["Z"][i];
     boundaries[i].a = geom["a"][i];
     boundaries[i].b = geom["b"][i];
     boundaries[i].c = geom["c"][i];
     boundaries[i].d = geom["d"][i];
     boundaries[i].plane_norm = geom["plane_norm"][i];
     boundaries[i].area = geom["area"][i];

     /*
     outfile << "geom(" << i+1 << ",:) = ["<<boundaries[i].x1 << ", " <<
        boundaries[i].z1 << ", " <<
        boundaries[i].x2 << ", " << boundaries[i].z2 << ", " <<
        boundaries[i].slope_dzdx << ", " << boundaries[i].intercept_z << ", " <<
        boundaries[i].length << ", " << boundaries[i].Z << "];" << std::endl;
     */
}   

  outfile.close();
  #else

  int nMaterials = geom["nMaterials"];
  std::cout << "nmat " << nMaterials << std::endl;
  for(int i=0 ; i<nLines ; i++)
  {
     boundaries[i].x1 = geom["x1"][i];
     boundaries[i].z1 = geom["z1"][i];
     boundaries[i].x2 = geom["x2"][i];
     boundaries[i].z2 = geom["z2"][i];
     std::cout << "z2 " << std::endl;
     boundaries[i].Z = geom["Z"][i];
     boundaries[i].slope_dzdx = geom["slope"][i];
     boundaries[i].intercept_z = geom["intercept"][i];
     boundaries[i].length = geom["length"][i];

     outfile << "geom(" << i+1 << ",:) = ["<<boundaries[i].x1 << ", " <<
        boundaries[i].z1 << ", " <<
        boundaries[i].x2 << ", " << boundaries[i].z2 << ", " <<
        boundaries[i].slope_dzdx << ", " << boundaries[i].intercept_z << ", " <<
        boundaries[i].length << ", " << boundaries[i].Z << "];" << std::endl;
  }   

  outfile.close();
    std::cout << "finished loop " << std::endl;
  boundaries[nLines].Z = geom["Z"][nLines];
  std::cout << " here 1" << std::endl;
  boundaries[nLines].y1 = geom["y1"];
  std::cout << " here 2" << std::endl;
  boundaries[nLines].y2 = geom["y2"];
  std::cout << " here 3" << std::endl;
  boundaries[nLines].periodic = geom["periodic"];
  std::cout << " here 4" << std::endl;
  #endif
  std::cout << "Starting Boundary Init..." << std::endl;

  //Applying background values at material boundaries
  std::for_each(boundaries.begin(), boundaries.end()-1,
            boundary_init(background_Z,background_amu,
            nR_Dens,nZ_Dens,DensGridr.data(),DensGridz.data(),ni.data(),
            nR_Bfield,nZ_Bfield,bfieldGridr.data(),
            bfieldGridz.data(),br.data(),bz.data(), bt.data(),
            nR_Temp,nZ_Temp,TempGridr.data(),
            TempGridz.data(),ti.data() ));

   std::cout << "Completed Boundary Init " << std::endl;
  //Efield
  #if USEPRESHEATHEFIELD > 0    
    #if PRESHEATH_INTERP == 0
      int nR_PreSheathEfield = 1;
      int nZ_PreSheathEfield = 1;
      sim::Array<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
      sim::Array<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), 
          PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
          PSEt(nR_PreSheathEfield*nZ_PreSheathEfield);
      PSEr[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Er");
      PSEz[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Ez");
      PSEt[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Et");
    #elif PRESHEATH_INTERP == 2
      int nR_PreSheathEfield;
      int nZ_PreSheathEfield;
      
      int e1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.gridNrString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.gridNzString"),
                  nR_PreSheathEfield,nZ_PreSheathEfield);
      
      sim::Array<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
      sim::Array<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), 
          PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
          PSEt(nR_PreSheathEfield*nZ_PreSheathEfield,0.0);
      
      int e2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.gridRString"), preSheathEGridr);
      
      int e3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.gridZString"), preSheathEGridz);
      
      int e4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.radialComponentString"), PSEr);
      
      int e5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.Efield.axialComponentString"), PSEz);
      
      //int e6 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
        //          cfg.lookup("backgroundPlasmaProfiles.Efield.toroidalComponentString"), PSEt);
    #endif
    
    std::string outnamePSEfieldR = "PSEfieldR.m";
    std::string outnamePSEfieldZ = "PSEfieldZ.m";
    std::string outnamePSEGridR = "PSEgridR.m";
    std::string outnamePSEGridZ = "PSEgridZ.m";
    OUTPUT1d(profiles_folder,outnamePSEGridR, nR_PreSheathEfield, &preSheathEGridr.front());
    OUTPUT1d(profiles_folder,outnamePSEGridZ, nZ_PreSheathEfield, &preSheathEGridz.front());
    OUTPUT2d(profiles_folder,outnamePSEfieldR, nR_PreSheathEfield, nZ_PreSheathEfield, &PSEr.front());
    OUTPUT2d(profiles_folder,outnamePSEfieldZ, nR_PreSheathEfield, nZ_PreSheathEfield, &PSEz.front());
  #else
    
      int nR_PreSheathEfield = 1;
      int nZ_PreSheathEfield = 1;
      sim::Array<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
      sim::Array<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), 
          PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
          PSEt(nR_PreSheathEfield*nZ_PreSheathEfield);
  #endif
    
  std::cout << "Completed presheath Efield Init " << std::endl;
  sim::Array<float> Efieldr(nR_Bfield*nZ_Bfield), Efieldz(nR_Bfield*nZ_Bfield),
                    Efieldt(nR_Bfield*nZ_Bfield),minDist(nR_Bfield*nZ_Bfield);

  #if USESHEATHEFIELD > 0
    #if EFIELD_INTERP == 1
      float thisE[3] = {0.0,0.0,0.0};
    
      for(int i=0;i<nR_Bfield;i++)
      {
         for(int j=0;j<nZ_Bfield;j++)
         {
             minDist[(nR_Bfield - 1 -i)*nZ_Bfield+(nZ_Bfield -1-j)] = 
                  getE ( bfieldGridr[i], 0.0, bfieldGridz[j],
                  thisE, boundaries.data(),nLines );
             Efieldr[i*nZ_Bfield+j] = thisE[0];
             Efieldz[i*nZ_Bfield+j] = thisE[2];
             Efieldt[i*nZ_Bfield+j] = thisE[1];
          }
      }
        
      int nR_closeGeom;
      int nZ_dtsEfield;
      
      int d1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNrString"),
                  cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNzString"),nR_dtsEfield,nZ_dtsEfield);
      
      sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
      sim::Array<float> dtsE(nR_dtsEfield*nZ_dtsEfield);
      
      int d2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridRString"), dtsEfieldGridr);
      
      std::cout << "got first grid " << dtsEfieldGridr.front() << std::endl;    
      int d3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridZString"), dtsEfieldGridz);
      
      std::cout << "got second grid" << dtsEfieldGridz.front() << std::endl;    
      
      int d4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                  cfg.lookup("backgroundPlasmaProfiles.dtsEfield.sheathDTS"), dtsE);
    #elif EFIELD_INTERP ==2
        int nR_dtsEfield, nZ_dtsEfield;
        
        int d1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                    cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNrString"),
                    cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridNzString"),
                    nR_dtsEfield,nZ_dtsEfield);
        
        sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
        sim::Array<float> dtsE(nR_dtsEfield*nZ_dtsEfield);
        
        int d2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                    cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridRString"), dtsEfieldGridr);
        
        int d3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                    cfg.lookup("backgroundPlasmaProfiles.dtsEfield.gridZString"), dtsEfieldGridz);
        
        int d4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.dtsEfield.fileString"),
                    cfg.lookup("backgroundPlasmaProfiles.dtsEfield.sheathDTS"), dtsE);
    #endif
  #else
    int nR_dtsEfield=1;
    int nZ_dtsEfield=1;
    sim::Array<float> dtsEfieldGridr(nR_dtsEfield), dtsEfieldGridz(nZ_dtsEfield);
    sim::Array<float> dtsE(nR_dtsEfield*nZ_dtsEfield);
  #endif

  std::string outnameEfieldR = "EfieldR.m";
  std::string outnameEfieldZ = "EfieldZ.m";
  std::string outnameEfieldT = "EfieldT.m";
  std::string outnameMinDist = "DistToSurface.m";
  OUTPUT2d(profiles_folder,outnameEfieldR, nR_Bfield, nZ_Bfield, &Efieldr.front());
  OUTPUT2d(profiles_folder,outnameEfieldZ, nR_Bfield, nZ_Bfield, &Efieldz.front());
  OUTPUT2d(profiles_folder,outnameEfieldT, nR_Bfield, nZ_Bfield, &Efieldt.front());
  OUTPUT2d(profiles_folder,outnameMinDist, nR_Bfield, nZ_Bfield, &minDist.front());


    int nR_closeGeom;
    int nZ_closeGeom;
    int n_closeGeomElements;
    int nY_closeGeom;
  #if GEOM_HASH > 0
#if USE3DTETGEOM >0
    int gi1 = read_profileNs(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridNrString"),
                cfg.lookup("geometry.gridNyString"),nR_closeGeom,nY_closeGeom);
    
    int gi2 = read_profileNs(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridNzString"),
                cfg.lookup("geometry.nearestNelementsString"),nZ_closeGeom,n_closeGeomElements);
std::cout << "3d tet geom hash " << nR_closeGeom << " " << nY_closeGeom << " " 
        << nZ_closeGeom << " " <<n_closeGeomElements << std::endl;
    sim::Array<float> closeGeomGridr(nR_closeGeom), closeGeomGridy(nY_closeGeom), closeGeomGridz(nZ_closeGeom);
    sim::Array<int> closeGeom(nR_closeGeom*nY_closeGeom*nZ_closeGeom*n_closeGeomElements);
    
    int gi3 = read_profile1d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridRString"), closeGeomGridr);
    
    int gi4 = read_profile1d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridYString"), closeGeomGridy);
    int gi5 = read_profile1d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridZString"), closeGeomGridz);
   
    int gi6 = read_profile3d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.closeGeomString"), closeGeom);
    std::cout << " finished geometry hashing read in " << std::endl;
#else
    int gi1 = read_profileNs(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridNrString"),
                cfg.lookup("geometry.gridNzString"),nR_closeGeom,nZ_closeGeom);
    
    int gi2 = read_profileNs(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridNrString"),
                cfg.lookup("geometry.nearestNelementsString"),nR_closeGeom,n_closeGeomElements);
    
    sim::Array<float> closeGeomGridr(nR_closeGeom),closeGeomGridy(1),
                         closeGeomGridz(nZ_closeGeom);
    sim::Array<int> closeGeom(nR_closeGeom*nZ_closeGeom*n_closeGeomElements);
    
    int gi3 = read_profile1d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridRString"), closeGeomGridr);
    
    int gi4 = read_profile1d(cfg.lookup("geometry.fileString"),
                cfg.lookup("geometry.gridZString"), closeGeomGridz);
   
    int gi5 = read_profile3d(cfg.lookup("geometry.fileString"),
             cfg.lookup("geometry.closeGeomString"), closeGeom);
#endif
#else
    nR_closeGeom = 1;
    nZ_closeGeom = 1;
    nY_closeGeom = 1;
    n_closeGeomElements = 1;
    sim::Array<float> closeGeomGridr(nR_closeGeom),closeGeomGridy(nY_closeGeom), closeGeomGridz(nZ_closeGeom);
    sim::Array<int> closeGeom(nR_closeGeom*nZ_closeGeom*n_closeGeomElements);
  #endif
    std::cout << "3d tet geom hash " << nR_closeGeom << " " << nY_closeGeom << " "
            << nZ_closeGeom << " " <<n_closeGeomElements << std::endl;
                
    int nR_closeGeom_sheath;
    int nZ_closeGeom_sheath;
    int n_closeGeomElements_sheath;
  
  #if GEOM_HASH_SHEATH  
    int gis1 = read_profileNs(cfg.lookup("geometry_sheath.fileString"),
                cfg.lookup("geometry_sheath.gridNrString"),
                cfg.lookup("geometry_sheath.gridNzString"),nR_closeGeom_sheath,nZ_closeGeom_sheath);
    
    int gis2 = read_profileNs(cfg.lookup("geometry_sheath.fileString"),
                cfg.lookup("geometry_sheath.gridNrString"),
                cfg.lookup("geometry_sheath.nearestNelementsString"),nR_closeGeom_sheath,n_closeGeomElements_sheath);
    
    sim::Array<float> closeGeomGridr_sheath(nR_closeGeom_sheath), 
                      closeGeomGridz_sheath(nZ_closeGeom_sheath);
    sim::Array<int> closeGeom_sheath(nR_closeGeom_sheath*nZ_closeGeom_sheath*n_closeGeomElements_sheath);
    
    int gis3 = read_profile1d(cfg.lookup("geometry_sheath.fileString"),
                cfg.lookup("geometry_sheath.gridRString"), closeGeomGridr_sheath);
    
    int gis4 = read_profile1d(cfg.lookup("geometry_sheath.fileString"),
                cfg.lookup("geometry_sheath.gridZString"), closeGeomGridz_sheath);
   
    int gis5 = read_profile3d(cfg.lookup("geometry_sheath.fileString"),
                cfg.lookup("geometry_sheath.closeGeomString"), closeGeom_sheath);
  #else
    nR_closeGeom_sheath = 1;
    nZ_closeGeom_sheath = 1;
    n_closeGeomElements_sheath = 1;
    sim::Array<float> closeGeomGridr_sheath(nR_closeGeom_sheath), 
                      closeGeomGridz_sheath(nZ_closeGeom_sheath);
    sim::Array<int> closeGeom_sheath(nR_closeGeom_sheath*nZ_closeGeom_sheath*n_closeGeomElements_sheath);
  #endif  

  #if SPECTROSCOPY > 0
    float netX0 = cfg.lookup("diagnostics.netx0");
    float netX1 = cfg.lookup("diagnostics.netx1");
    float netY0 = cfg.lookup("diagnostics.nety0");
    float netY1 = cfg.lookup("diagnostics.nety1");
    float netZ0 = cfg.lookup("diagnostics.netz0");
    float netZ1 = cfg.lookup("diagnostics.netz1");
    int net_nX = cfg.lookup("diagnostics.nX");
    int net_nY = cfg.lookup("diagnostics.nY");
    int net_nZ = cfg.lookup("diagnostics.nZ");
    Setting& diagn = cfg.lookup("diagnostics");
    int nBins = cfg.lookup("diagnostics.densityChargeBins");//1;//diagn["densityChargeBins"].getLength();
  
    #if USECYLSYMM > 0
      sim::Array<float> net_Bins((nBins+1)*net_nX*net_nZ);
      /*
      for (int i=0; i<nBins*net_nX*net_nZ; i++)
          {
              std::cout << "i " << i << std::endl;
            net_Bins[i] = 0;
              std::cout << "net bins " << net_Bins[i] << std::endl;
            
          }
      */
      sim::Array<float> gridX_bins(net_nX), gridZ_bins(net_nZ);

      for (int i=0; i< net_nX ; i++)
      {
         gridX_bins[i] = netX0 + 1.0/(net_nX-1)*i*(netX1-netX0);
      }

      for (int i=0; i< net_nZ ; i++)
      {
         gridZ_bins[i] = netZ0 + i*1.0/(net_nZ-1)*(netZ1-netZ0);
      }
    #endif
  #endif    

  // Perp DiffusionCoeff initialization - only used when Diffusion interpolator is = 0
  float perpDiffusionCoeff=0.0;
  if(world_rank ==0)
  {
  if (cfg.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",perpDiffusionCoeff))
  {}
  else
  {std::cout << "ERROR: could not get perpendicular diffusion coefficient from input file" << std::endl;}
  }
  #if USE_MPI > 0
      MPI_Bcast(&perpDiffusionCoeff,1,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
  #endif
//Surface model import
  int nE_sputtRefCoeff = 1, nA_sputtRefCoeff=1;
  int nE_sputtRefDistIn = 1, nA_sputtRefDistIn=1;
  int nE_sputtRefDistOut = 1, nA_sputtRefDistOut = 1;
  int nE_sputtRefDistOutRef = 1,nDistE_surfaceModelRef = 1;
  int nDistE_surfaceModel = 1,nDistA_surfaceModel = 1;
#if USESURFACEMODEL > 0
  std::string surfaceModelCfg = "surfaceModel.";
  std::string surfaceModelFile;
  #if USE_MPI > 0 
    if(world_rank == 0)
    {
  #endif
  getVariable(cfg,surfaceModelCfg+"fileString",surfaceModelFile);
  nE_sputtRefCoeff = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefCoeffString");
  nA_sputtRefCoeff = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefCoeffString");
  nE_sputtRefDistIn = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefDistInString");
  nA_sputtRefDistIn = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefDistInString");
  nE_sputtRefDistOut = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefDistOutString");
  nE_sputtRefDistOutRef = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefDistOutStringRef");
  nA_sputtRefDistOut = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefDistOutString");
  nDistE_surfaceModel = nE_sputtRefDistIn*nA_sputtRefDistIn*nE_sputtRefDistOut;
  nDistE_surfaceModelRef = nE_sputtRefDistIn*nA_sputtRefDistIn*nE_sputtRefDistOutRef;
  nDistA_surfaceModel = nE_sputtRefDistIn*nA_sputtRefDistIn*nA_sputtRefDistOut;
  std::cout <<  " got dimensions of surface model " << std::endl;
  #if USE_MPI > 0
    }
      MPI_Bcast(&nE_sputtRefCoeff, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nA_sputtRefCoeff, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nE_sputtRefDistIn, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nA_sputtRefDistIn, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nE_sputtRefDistOut, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nE_sputtRefDistOutRef, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nA_sputtRefDistOut, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nDistE_surfaceModel, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nDistE_surfaceModelRef, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nDistA_surfaceModel, 1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
  #endif
#endif
  sim::Array<float> E_sputtRefCoeff(nE_sputtRefCoeff), A_sputtRefCoeff(nA_sputtRefCoeff),
                    Elog_sputtRefCoeff(nE_sputtRefCoeff),
                    energyDistGrid01(nE_sputtRefDistOut),
                    energyDistGrid01Ref(nE_sputtRefDistOutRef),
                    angleDistGrid01(nA_sputtRefDistOut),
                    spyl_surfaceModel(nE_sputtRefCoeff*nA_sputtRefCoeff),
                    rfyl_surfaceModel(nE_sputtRefCoeff*nA_sputtRefCoeff),
                    E_sputtRefDistIn(nE_sputtRefDistIn), A_sputtRefDistIn(nA_sputtRefDistIn),
                    Elog_sputtRefDistIn(nE_sputtRefDistIn),
                    E_sputtRefDistOut(nE_sputtRefDistOut), 
                    E_sputtRefDistOutRef(nE_sputtRefDistOutRef),
		    Aphi_sputtRefDistOut(nA_sputtRefDistOut),Atheta_sputtRefDistOut(nA_sputtRefDistOut),
                    AphiDist_Y(nDistA_surfaceModel),AthetaDist_Y(nDistA_surfaceModel),
		    EDist_Y(nDistE_surfaceModel),
                    AphiDist_R(nDistA_surfaceModel),AthetaDist_R(nDistA_surfaceModel),
		    EDist_R(nDistE_surfaceModelRef),
                    AphiDist_CDF_Y(nDistA_surfaceModel),AthetaDist_CDF_Y(nDistA_surfaceModel),
		    EDist_CDF_Y(nDistE_surfaceModel),
                    AphiDist_CDF_R(nDistA_surfaceModel),AthetaDist_CDF_R(nDistA_surfaceModel),
		    EDist_CDF_R(nDistE_surfaceModelRef),
                    AphiDist_CDF_Y_regrid(nDistA_surfaceModel),AthetaDist_CDF_Y_regrid(nDistA_surfaceModel),
		    EDist_CDF_Y_regrid(nDistE_surfaceModel),
                    AphiDist_CDF_R_regrid(nDistA_surfaceModel),AthetaDist_CDF_R_regrid(nDistA_surfaceModel),
		    EDist_CDF_R_regrid(nDistE_surfaceModelRef);
#if USESURFACEMODEL > 0
  #if USE_MPI > 0 
    if(world_rank == 0)
    {
  #endif
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefCoeff",E_sputtRefCoeff[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"A_sputtRefCoeff",A_sputtRefCoeff[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefDistIn",E_sputtRefDistIn[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"A_sputtRefDistIn",A_sputtRefDistIn[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefDistOut",E_sputtRefDistOut[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefDistOutRef",E_sputtRefDistOutRef[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"Aphi_sputtRefDistOut",Aphi_sputtRefDistOut[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"Atheta_sputtRefDistOut",Atheta_sputtRefDistOut[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"sputtYldString",spyl_surfaceModel[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"reflYldString",rfyl_surfaceModel[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"EDist_Y",EDist_Y[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"AphiDist_Y",AphiDist_Y[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"AthetaDist_Y",AthetaDist_Y[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"EDist_R",EDist_R[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"AphiDist_R",AphiDist_R[0]);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"AthetaDist_R",AthetaDist_R[0]);
  //for(int i=0;i<nDistE_surfaceModel;i++)
  //{
  //    std::cout << " Edist diff Y " << EDist_Y[i] << " " << EDist_R[i] << std::endl;
  //}
  for(int i=0;i<nE_sputtRefCoeff;i++)
  {
      Elog_sputtRefCoeff[i] = log10(E_sputtRefCoeff[i]);
      std::cout << " EsputtRefCoeff and Elog " << E_sputtRefCoeff[i] << " " << Elog_sputtRefCoeff[i] << std::endl;
  }
  for(int i=0;i<nE_sputtRefDistIn;i++)
  {
      Elog_sputtRefDistIn[i] = log10(E_sputtRefDistIn[i]);
  }
  for(int i=0;i<nE_sputtRefDistOut;i++)
  {
      energyDistGrid01[i] = i*1.0/nE_sputtRefDistOut;
  }
  for(int i=0;i<nE_sputtRefDistOutRef;i++)
  {
      energyDistGrid01Ref[i] = i*1.0/nE_sputtRefDistOutRef;
  }
  for(int i=0;i<nA_sputtRefDistOut;i++)
  {
     angleDistGrid01[i] = i*1.0/nA_sputtRefDistOut;
     //std::cout << " angleDistGrid01[i] " << angleDistGrid01[i] << std::endl;
  }
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,EDist_Y.data(),EDist_CDF_Y.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,AphiDist_Y.data(),AphiDist_CDF_Y.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,AthetaDist_Y.data(),AthetaDist_CDF_Y.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOutRef,EDist_R.data(),EDist_CDF_R.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,AphiDist_R.data(),AphiDist_CDF_R.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,AthetaDist_R.data(),AthetaDist_CDF_R.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,AthetaDist_R.data(),AthetaDist_CDF_R.data());
  //for(int k=0;k<nE_sputtRefDistOut;k++)
  //{
  //      std::cout << "Edist_CDF_Y " << EDist_CDF_Y[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
  ////      std::cout << "cosDist_CDFR " << EDist_CDF_R[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
  //}
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,Aphi_sputtRefDistOut[nA_sputtRefDistOut-1],AphiDist_CDF_Y.data(),AphiDist_CDF_Y_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,Atheta_sputtRefDistOut[nA_sputtRefDistOut-1],AthetaDist_CDF_Y.data(),AthetaDist_CDF_Y_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,energyDistGrid01.data(),nE_sputtRefDistOut,E_sputtRefDistOut[nE_sputtRefDistOut-1],EDist_CDF_Y.data(),EDist_CDF_Y_regrid.data());
 //std::cout << "max value " << E_sputtRefDistOut[nE_sputtRefDistOut-1] << std::endl;
 //for(int k=0;k<60;k++)
 // {
 //     std::cout << "Edis amd cdf " << k << " " << EDist_CDF_Y[k] << " " <<EDist_CDF_Y_regrid[k] << std::endl; 
 // }
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,Aphi_sputtRefDistOut[nA_sputtRefDistOut-1],AphiDist_CDF_R.data(),AphiDist_CDF_R_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,Atheta_sputtRefDistOut[nA_sputtRefDistOut-1],AthetaDist_CDF_R.data(),AthetaDist_CDF_R_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOutRef,energyDistGrid01Ref.data(),nE_sputtRefDistOutRef,E_sputtRefDistOutRef[nE_sputtRefDistOutRef-1],EDist_CDF_R.data(),EDist_CDF_R_regrid.data());
 // regrid2dCDF(nE_surfaceModel,nA_surfaceModel,nEdistBins_surfaceModel,energyDistGrid01.data(),nEdistBins_surfaceModel,100.0,energyDist_CDF.data(),energyDist_CDFregrid.data());
 //for(int k=0;k<nE_sputtRefDistOut;k++)
 // {
 //       std::cout << "EDist_CDFregridY " << EDist_CDF_Y_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
 ////       std::cout << "cosDist_CDFregridR " << EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
 // }
 //for(int k=0;k<nA_sputtRefDistOut;k++)
 // {
 //       std::cout << "ADist_CDFregridY " << k << " " << AphiDist_Y[k]<< " " << AphiDist_CDF_Y[0*nA_sputtRefDistIn*nA_sputtRefDistOut + 0*nA_sputtRefDistOut+k]<< " " <<  AphiDist_CDF_Y_regrid[0*nA_sputtRefDistIn*nA_sputtRefDistOut + 0*nA_sputtRefDistOut+k] << std::endl;
 ////       std::cout << "cosDist_CDFregridR " << EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
 // }
 //for(int k=0;k<nA_sputtRefDistOut;k++)
 // {
 //       std::cout << "ADist_CDFregridR " << k << " " << AthetaDist_R[k]<< " " << AthetaDist_CDF_R[0*nA_sputtRefDistIn*nA_sputtRefDistOut + 0*nA_sputtRefDistOut+k]<< " " <<  AthetaDist_CDF_R_regrid[0*nA_sputtRefDistIn*nA_sputtRefDistOut + 0*nA_sputtRefDistOut+k] << std::endl;
 ////       std::cout << "cosDist_CDFregridR " << EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
 // }
  //float spylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff, nE_sputtRefCoeff,A_sputtRefCoeff.data(),
    //                          Elog_sputtRefCoeff.data(),spyl_surfaceModel.data());
  //float rfylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff, nE_sputtRefCoeff,A_sputtRefCoeff.data(),
      //                        Elog_sputtRefCoeff.data(),rfyl_surfaceModel.data());
  float spylAInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,AphiDist_CDF_Y_regrid.data() );
  float spylAthetaInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,AthetaDist_CDF_Y_regrid.data() );
 float sputEInterpVal = interp3d ( 0.44,63.0,log10(10.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
              energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,EDist_CDF_Y_regrid.data() );
  float rfylAInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,AphiDist_CDF_R_regrid.data() );
  float rfylAthetaInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,AthetaDist_CDF_R_regrid.data() );
 float rflEInterpVal = interp3d ( 0.44,63.0,log10(10.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
              energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,EDist_CDF_R_regrid.data() );
  //float rflAInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
   //       angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,ADist_CDF_R_regrid.data() );
 //float rflEInterpVal = interp3d ( 0.44,5.0,log10(250.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
     //         energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,EDist_CDF_R_regrid.data() );
  //std::cout << "Finished surface model import " <<spylInterpVal << " " <<  spylAInterpVal << " " << sputEInterpVal << " "<< rfylInterpVal<< " " << rflAInterpVal << " " << rflEInterpVal <<  std::endl; 
  std::cout << "Finished surface model import sputtering"  << spylAInterpVal << " " << spylAthetaInterpVal << " " << sputEInterpVal  <<  std::endl; 
  std::cout << "Finished surface model import reflection"  << rfylAInterpVal << " " << rfylAthetaInterpVal << " " << rflEInterpVal  <<  std::endl; 
  #if USE_MPI > 0
    }
      MPI_Bcast(E_sputtRefCoeff.data(), nE_sputtRefCoeff,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(A_sputtRefCoeff.data(), nA_sputtRefCoeff,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(Elog_sputtRefCoeff.data(), nE_sputtRefCoeff,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(energyDistGrid01.data(), nE_sputtRefDistOut,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(energyDistGrid01Ref.data(), nE_sputtRefDistOutRef,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(angleDistGrid01.data(), nA_sputtRefDistOut,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(spyl_surfaceModel.data(), nE_sputtRefCoeff*nA_sputtRefCoeff,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(rfyl_surfaceModel.data(), nE_sputtRefCoeff*nA_sputtRefCoeff,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(E_sputtRefDistIn.data(), nE_sputtRefDistIn,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(A_sputtRefDistIn.data(), nA_sputtRefDistIn,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(Elog_sputtRefDistIn.data(), nE_sputtRefDistIn,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(E_sputtRefDistOut.data(), nE_sputtRefDistOut,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(E_sputtRefDistOutRef.data(), nE_sputtRefDistOutRef,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(Aphi_sputtRefDistOut.data(), nA_sputtRefDistOut,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(Atheta_sputtRefDistOut.data(), nA_sputtRefDistOut,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_Y.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_Y.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_Y.data(), nDistE_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_R.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_R.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_R.data(), nDistE_surfaceModelRef,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_CDF_Y.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_CDF_Y.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_CDF_Y.data(), nDistE_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_CDF_R.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_CDF_R.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_CDF_R.data(), nDistE_surfaceModelRef,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_CDF_Y_regrid.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_CDF_Y_regrid.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_CDF_Y_regrid.data(), nDistE_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AphiDist_CDF_R_regrid.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(AthetaDist_CDF_R_regrid.data(), nDistA_surfaceModel,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(EDist_CDF_R_regrid.data(), nDistE_surfaceModelRef,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
  #endif
#endif
  // Particle time stepping control
  int ionization_nDtPerApply  = cfg.lookup("timeStep.ionization_nDtPerApply");
  int collision_nDtPerApply  = cfg.lookup("timeStep.collision_nDtPerApply");

  #ifdef __CUDACC__
    cout<<"Using THRUST"<<endl;
  #else
    cout<<"Not using THRUST"<<endl;
  #endif

  float dt = cfg.lookup("timeStep.dt");
  //float nPtsPerGyroOrbit = cfg.lookup("timeStep.nPtsPerGyroOrbit");
  //dt = 2.4e-7/100.0;

  const int nP = cfg.lookup("impurityParticleSource.nP");
  cout << "Number of particles: " << nP << endl;              
  long nParticles = nP;
  int nT = cfg.lookup("timeStep.nT");
  cout << "Number of time steps: " << nT << " With dt = " << dt << endl; 

#if PARTICLE_SOURCE == 0
float x = cfg.lookup("impurityParticleSource.initialConditions.x_start");
float y = cfg.lookup("impurityParticleSource.initialConditions.y_start");
float z = cfg.lookup("impurityParticleSource.initialConditions.z_start");

float Ex = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
float Ey = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
float Ez = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");

float amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
float Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");
float charge = cfg.lookup("impurityParticleSource.initialConditions.charge");
//    Particle p1(x,y,z,Ex,Ey,Ez,Z,amu,charge);
//    sim::Array<Particle> particleArray(nParticles,p1);
auto particleArray = new Particles(nParticles);
for (int i=0; i< nP ; i++)
{
particleArray->setParticle(i,x, y, z, Ex, Ey, Ez, Z, amu, charge);
}
#elif PARTICLE_SOURCE == 1
float x;
float y;
float z;

float Ex;
float Ey;
float Ez;

float amu;
float Z;
float charge;
    float impurity_Z = cfg.lookup("impurityParticleSource.Z");
    int nImpurityBoundaries = 0;
    for (int i=0; i<nLines;i++)
    {
        if(boundaries[i].Z == impurity_Z)
        {
            nImpurityBoundaries++;
        }
    }
    std::cout << "n Impurity Boundaries to launch from " << nImpurityBoundaries << std::endl;
    std::vector<int> boundaryIndex_ImpurityLaunch(nImpurityBoundaries);

    int count = 0;
    for (int i=0; i<nLines;i++)
    {
        if(boundaries[i].Z == impurity_Z)
        {
            boundaryIndex_ImpurityLaunch[count] = i;
            count++;
            std::cout << "Boundary indices " << i << std::endl;
        }
    }
    
    int impuritiesPerBoundary = nP/nImpurityBoundaries;
      
    //sim::Array<Particle> particleArray(nParticles);  
auto particleArray = new Particles(nParticles);

      std::uniform_real_distribution<float> distributionForSeeds(0,1e6);
#if FIXEDSEEDS ==0
    std::random_device randDevice;
    std::default_random_engine generator0(randDevice());
#else
    float randDevice = 6.5298E+5;
    std::default_random_engine generator0(randDevice);
#endif
    
    std::vector<float> boundarySeeds0(4*nImpurityBoundaries);
    std::generate( boundarySeeds0.begin(), boundarySeeds0.end(), [&]() { return distributionForSeeds(generator0); } );
    std::uniform_real_distribution<float> dist01(0.0, 1.0);
    float rand0 = 0.0;
    float rand1 = 0.0;
    float rand2 = 0.0;
    float rand3 = 0.0;

    std::vector<std::mt19937> s0(4*nImpurityBoundaries);
    
    float E0 = 0.0;
//Create Thompson Distribution
    float surfaceBindingEnergy = cfg.lookup("impurityParticleSource.source_material_SurfaceBindingEnergy");
    std::cout << "surface binding energy " << surfaceBindingEnergy << std::endl;
    int nThompDistPoints = 200;
    float max_Energy = 100.0;
    std::vector<float> ThompsonDist(nThompDistPoints),CumulativeDFThompson(nThompDistPoints);
    for(int i=0;i<nThompDistPoints;i++)
        {
            ThompsonDist[i] = (i*max_Energy/nThompDistPoints)/pow((i*max_Energy/nThompDistPoints) + surfaceBindingEnergy,3);
            if(i==0)
            {
                CumulativeDFThompson[i] = ThompsonDist[i]; 
            }
            else
            {
                CumulativeDFThompson[i] = CumulativeDFThompson[i-1]+ThompsonDist[i];
            }
        }
    for(int i=0;i<nThompDistPoints;i++)
        {
            CumulativeDFThompson[i] = CumulativeDFThompson[i]/CumulativeDFThompson[nThompDistPoints-1];
            //std::cout << "energy and CDF" << i*max_Energy/nThompDistPoints << " " << CumulativeDFThompson[i] << std::endl;
        }

    for(int j=0; j<4*nImpurityBoundaries;j++)
        {
            std::mt19937  s(boundarySeeds0[j]);
            s0[j] = s;
        }
    // Particle p1(0.0,0.0,0.0,0.0,0.0,0.0,0,0.0);
    for (int i=0; i< nImpurityBoundaries;i++)
    {
        for(int j=0; j<impuritiesPerBoundary; j++)
        {
            //Set boundary interval, properties, and random number gen
        if (i==0)
        {
            rand0 = dist01(s0[0]);
            x = boundaries[boundaryIndex_ImpurityLaunch[i]].x1 + 
                boundaries[boundaryIndex_ImpurityLaunch[i]].length*rand0;//1.4290;
            //std::cout << "start pos 1 " << x << std::endl;
            z = -1.2540+0.00001;
            rand1 = dist01(s0[1]);
            rand2 = dist01(s0[2]);
            rand3 = dist01(s0[3]);
            E0 = interp1dUnstructured(rand2,nThompDistPoints, max_Energy, &CumulativeDFThompson.front());
            Ex = E0*cos(3.1415*rand1)*sin(3.1415*rand3);
            Ey = E0*cos(3.1415*rand3);
            Ez = E0*sin(3.1415*rand1)*sin(3.1415*rand3);
        }
        else
        {
            rand0 = dist01(s0[4]);
            x = boundaries[boundaryIndex_ImpurityLaunch[i]].x1 + boundaries[boundaryIndex_ImpurityLaunch[i]].length*rand0;
            //x = 1.3450;
            //std::cout << "start pos 2 " << x << std::endl;
            z = -1.3660+0.00001;
            rand1 = dist01(s0[5]);
            rand2 = dist01(s0[6]);
            rand3 = dist01(s0[7]);
            E0 = interp1dUnstructured(rand2,nThompDistPoints, max_Energy, &CumulativeDFThompson.front());
            Ex = E0*cos(3.1415*rand1)*sin(3.1415*rand3);
            Ey = E0*cos(3.1415*rand3);
            Ez = E0*sin(3.1415*rand1)*sin(3.1415*rand3);
        }
particleArray->setParticle((i * impuritiesPerBoundary + j),x, 0.0, z, Ex, Ey, Ez, 74, 184.0, charge);            
        //Particle p1(x,0.0,z,Ex,Ey,Ez,74,184.0,charge);
            //particleArray[i*impuritiesPerBoundary + j] = p1;
          //  particleArray[i*impuritiesPerBoundary + j] = p1;
            //std::cout << " E0 " << E0 << std::endl;
            //std::cout << "vy " << particleArray[i*impuritiesPerBoundary + j].vy << " " << Ey << std::endl;
            //std::cout << "vx " << particleArray[i*impuritiesPerBoundary + j].vx << " " << Ex << std::endl;
            //std::cout << "vz " << particleArray[i*impuritiesPerBoundary + j].vz << " " << Ez << std::endl;
        }
    }
#endif


#if GEOM_TRACE > 0       
    std::uniform_real_distribution<float> dist2(0,1);
    //std::random_device rd2;
    //std::default_random_engine generator2(rd2());
        float randDevice02 = 6.52E+5;
        std::default_random_engine generator2(randDevice02);
    std::cout << "Randomizing velocities to trace geometry. " << std::endl;

    for (int i=0 ; i<nParticles ; i++)
    {   float theta = dist2(generator2)*2*3.1415;
        float phi = dist2(generator2)*3.1415;
        float mag = 2e3;
        particleArray->vx[i] = mag*cos(theta)*sin(phi);
        particleArray->vy[i] = mag*sin(theta)*sin(phi);
        particleArray->vz[i] = mag*cos(phi);
    }
#endif

#if PARTICLE_TRACKS > 0
    int subSampleFac = 1;
    float subSampleFacf = 1.0;
#if USE_CUDA > 0
    sim::Array<float> positionHistoryX(nP*nT/subSampleFac);
    sim::Array<float> positionHistoryY(nP*nT/subSampleFac);
    sim::Array<float> positionHistoryZ(nP*nT/subSampleFac);
    sim::Array<float> velocityHistoryX(nP*nT/subSampleFac);
    sim::Array<float> velocityHistoryY(nP*nT/subSampleFac);
    sim::Array<float> velocityHistoryZ(nP*nT/subSampleFac);
    sim::Array<float> chargeHistory(nP*nT/subSampleFac);
#else
    float **positionHistoryX;
    float **positionHistoryY;
    float **positionHistoryZ;
    float **velocityHistoryX;
    float **velocityHistoryY;
    float **velocityHistoryZ;
    float **chargeHistory;
    positionHistoryX = new float* [nP];
    positionHistoryY = new float* [nP];
    positionHistoryZ = new float* [nP];
    velocityHistoryX = new float* [nP];
    velocityHistoryY = new float* [nP];
    velocityHistoryZ = new float* [nP];
    chargeHistory = new float* [nP];
    positionHistoryX[0] = new float [nT*nP/subSampleFac];
    positionHistoryY[0] = new float [nT*nP/subSampleFac];
    positionHistoryZ[0] = new float [nT*nP/subSampleFac];
    velocityHistoryX[0] = new float [nT*nP/subSampleFac];
    velocityHistoryY[0] = new float [nT*nP/subSampleFac];
    velocityHistoryZ[0] = new float [nT*nP/subSampleFac];
    chargeHistory[0] = new float [nT*nP/subSampleFac];
    for(int i=0 ; i<nP ; i++)
    {
        positionHistoryX[i] = &positionHistoryX[0][i*nT/subSampleFac];
        positionHistoryY[i] = &positionHistoryY[0][i*nT/subSampleFac];
        positionHistoryZ[i] = &positionHistoryZ[0][i*nT/subSampleFac];
        velocityHistoryX[i] = &velocityHistoryX[0][i*nT/subSampleFac];
        velocityHistoryY[i] = &velocityHistoryY[0][i*nT/subSampleFac];
        velocityHistoryZ[i] = &velocityHistoryZ[0][i*nT/subSampleFac];
        chargeHistory[i] = &chargeHistory[0][i*nT/subSampleFac];
        for(int j=0 ; j<nT/subSampleFac ; j++)
        {
            positionHistoryX[i][j] = 0.0;
            positionHistoryY[i][j] = 0.0;
            positionHistoryZ[i][j] = 0.0;
            velocityHistoryX[i][j] = 0.0;
            velocityHistoryY[i][j] = 0.0;
            velocityHistoryZ[i][j] = 0.0;
            chargeHistory[i][j] = 0.0;
        }
    }
    const int* displ=&pDisplacement[0];
    const int* phpn=&pHistPerNode[0];
    std::cout << "history array length " << nHistories << std::endl;
    #if USE_CUDA > 0
      sim::Array<float> positionHistoryX(nHistories);
      sim::Array<float> positionHistoryXgather(nHistories,0.0);
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
      std::vector<float> positionHistoryXgather(nHistories,0.0);
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
  float* finalPosX = new float[nP];
  float* finalPosY = new float[nP];
  float* finalPosZ = new float[nP];
  float* finalVx = new float[nP];
  float* finalVy = new float[nP];
  float* finalVz = new float[nP];
  float* transitTime = new float[nP];
  float* hitWall = new float[nP];
  #if USE_BOOST
    //cpu_timer timer;
  #endif

std::uniform_real_distribution<float> dist(0,1e6);

#if FIXEDSEEDS == 0
    std::random_device rd;
    std::default_random_engine generator(rd());
#endif
  thrust::counting_iterator<std::size_t> particleBegin(0);  
  thrust::counting_iterator<std::size_t> particleEnd(nParticles-1);
  thrust::counting_iterator<std::size_t> particleOne(1);
    auto randInitStart_clock = Time::now();
    
  #if PARTICLESEEDS > 0
    #if USE_CUDA  
      sim::Array<curandState> state1(nParticles);
    #else
      sim::Array<std::mt19937> state1(nParticles);
    #endif
    #if USEIONIZATION > 0 || USERECOMBINATION > 0 || USEPERPDIFFUSION > 0 || USECOULOMBCOLLISIONS > 0 || USESURFACEMODEL > 0
#if USE_CUDA
      thrust::for_each(thrust::device, particleBegin,particleEnd,
      curandInitialize(&state1[0],0));
#else
      std::random_device randDeviceInit; 
      std::cout << "Initializing curand seeds " << std::endl;
      //thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,
      //                     curandInitialize(&state1[0],randDeviceInit,0));
       //std::mt19937 s0(randDeviceInit);
       for(int i=world_rank*nP/world_size;i<(world_rank+1)*nP/world_size;i++)
       {
           std::mt19937 s0(randDeviceInit());
           state1[i] = s0;
       }
#endif
      #if USE_CUDA
        cudaDeviceSynchronize();
      #endif
    #endif
  #endif
    auto randInitEnd_clock = Time::now();
    fsec fsRandInit = randInitEnd_clock - randInitStart_clock;
    printf("Random Number Initialize time for node %i          is %6.3f (secs) \n", world_rank,fsRandInit.count());

    float moveTime = 0.0;
    float geomCheckTime = 0.0;
    float ionizTime = 0.0;
    int* dev_tt;
    cudaMallocManaged(&dev_tt, sizeof(int));
    int tt=0;
    move_boris move_boris0(particleArray,dt,boundaries.data(), nLines,
        nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
        &br.front(),&bz.front(),&by.front(),
        nR_PreSheathEfield,nY_PreSheathEfield,nZ_PreSheathEfield,
        &preSheathEGridr.front(),&preSheathEGridy.front(),&preSheathEGridz.front(),
        &PSEr.front(),&PSEz.front(),&PSEt.front(),
            nR_closeGeom_sheath,nY_closeGeom_sheath,nZ_closeGeom_sheath,n_closeGeomElements_sheath,
            &closeGeomGridr_sheath.front(),&closeGeomGridy_sheath.front(),&closeGeomGridz_sheath.front(),
            &closeGeom_sheath.front());
     geometry_check geometry_check0(particleArray,nLines,&boundaries[0],surfaces,dt,
                        nHashes,nR_closeGeom.data(),nY_closeGeom.data(),nZ_closeGeom.data(),n_closeGeomElements.data(),
                        &closeGeomGridr.front(),&closeGeomGridy.front(),&closeGeomGridz.front(),
                        &closeGeom.front(),
                        nEdist, E0dist, Edist, nAdist, A0dist, Adist);
      #if USE_SORT > 0
        sortParticles sort0(particleArray,nP,0.001,dev_tt,10000,pStartIndx.data(),nActiveParticlesOnRank.data(),world_rank,&state1.front());
      #endif
#if SPECTROSCOPY > 0
            thrust::for_each(thrust::device, particleBegin,particleEnd,
                    spec_bin(particleArray,nBins,net_nX, net_nZ, &gridX_bins.front(),
                        &gridZ_bins.front(), &net_Bins.front(),dt) );
#endif            
#if USEIONIZATION > 0
        thrust::for_each(thrust::device, particleBegin,particleEnd,
                ionize(particleArray, dt,&state1.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),  
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                    nTemperaturesIonize, nDensitiesIonize,&gridTemperature_Ionization.front(),
                    &gridDensity_Ionization.front(), &rateCoeff_Ionization.front());
#endif
#if USERECOMBINATION > 0
                recombine recombine0(particleArray, dt,&state1.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),  
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                    nTemperaturesRecombine,nDensitiesRecombine,
                    gridTemperature_Recombination.data(),gridDensity_Recombination.data(),
                    rateCoeff_Recombination.data());
#endif
#if USEPERPDIFFUSION > 0
        thrust::for_each(particleArray.begin(), particleArray.end(),
                crossFieldDiffusion(dt,perpDiffusionCoeff,
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));
            
        thrust::for_each(particleArray.begin(), particleArray.end(),
                    geometry_check(nLines,BoundaryDevicePointer,dt,tt) );
#endif
#if USECOULOMBCOLLISIONS > 0
        thrust::for_each(particleArray.begin(), particleArray.end(), 
                coulombCollisions(dt,
                    nR_flowV,nZ_flowV,&flowVGridr.front(),&flowVGridz.front(),
                    &flowVr.front(),&flowVz.front(),&flowVt.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),    
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front()
                    background_Z,background_amu, 
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));

#endif
#if USETHERMALFORCE > 0
        thrust::for_each(particleArray.begin(), particleArray.end(),
                thermalForce(dt,background_amu,
                    nR_gradT,nZ_gradT,GradTGridRDevicePointer,GradTGridZDevicePointer,
                    GradTiRDevicePointer,GradTiZDevicePointer, GradTiTDevicePointer, 
                    GradTeRDevicePointer, GradTeZDevicePointer, GradTeTDevicePointer, 
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));
#endif

#if USESURFACEMODEL > 0
                reflection reflection0(particleArray,dt,&state1.front(),nLines,&boundaries[0],surfaces,
                    nE_sputtRefCoeff, nA_sputtRefCoeff,A_sputtRefCoeff.data(),
                    Elog_sputtRefCoeff.data(),spyl_surfaceModel.data(), rfyl_surfaceModel.data(),
                    nE_sputtRefDistOut,nE_sputtRefDistOutRef, nA_sputtRefDistOut,nE_sputtRefDistIn,nA_sputtRefDistIn,
                    Elog_sputtRefDistIn.data(),A_sputtRefDistIn.data(),
                    E_sputtRefDistOut.data(),E_sputtRefDistOutRef.data(),Aphi_sputtRefDistOut.data(),
                    energyDistGrid01.data(),energyDistGrid01Ref.data(),angleDistGrid01.data(),
                    EDist_CDF_Y_regrid.data(),AphiDist_CDF_Y_regrid.data(),
                    EDist_CDF_R_regrid.data(),AphiDist_CDF_R_regrid.data(),
                    nEdist, E0dist, Edist, nAdist, A0dist, Adist) ;
#endif        

#if PARTICLE_TRACKS >0
      history history0(particleArray,tt,nT,subSampleFac,nP,&positionHistoryX.front(),
      &positionHistoryY.front(),&positionHistoryZ.front(),&velocityHistory.front(),
      &velocityHistoryX.front(),&velocityHistoryY.front(),
      &velocityHistoryZ.front(),&chargeHistory.front(),
      &weightHistory.front());
#endif
  #if FORCE_EVAL > 0
    if(world_rank ==0)
    {
      int nR_force, nZ_force;
      float forceX0, forceX1,forceZ0,forceZ1,testEnergy;
      std::string forceCfg = "forceEvaluation.";
      
      getVariable(cfg,forceCfg+"nR",nR_force);
      getVariable(cfg,forceCfg+"nZ",nZ_force);
      std::vector<float> forceR(nR_force,0.0),forceZ(nZ_force,0.0);
      std::vector<float> dvEr(nR_force*nZ_force,0.0),dvEz(nR_force*nZ_force,0.0),dvEt(nR_force*nZ_force,0.0);
      std::vector<float> dvBr(nR_force*nZ_force,0.0),dvBz(nR_force*nZ_force,0.0),dvBt(nR_force*nZ_force,0.0);
      std::vector<float> dvCollr(nR_force*nZ_force,0.0),dvCollz(nR_force*nZ_force,0.0),dvCollt(nR_force*nZ_force,0.0);
      std::vector<float> dvITGr(nR_force*nZ_force,0.0),dvITGz(nR_force*nZ_force,0.0),dvITGt(nR_force*nZ_force,0.0);
      std::vector<float> dvETGr(nR_force*nZ_force,0.0),dvETGz(nR_force*nZ_force,0.0),dvETGt(nR_force*nZ_force,0.0);
      getVariable(cfg,forceCfg+"X0",forceX0);
      getVariable(cfg,forceCfg+"X1",forceX1);
      getVariable(cfg,forceCfg+"Z0",forceZ0);
      getVariable(cfg,forceCfg+"Z1",forceZ1);
      getVariable(cfg,forceCfg+"particleEnergy",testEnergy);
      for(int i=0;i<nR_force;i++)
      {
        forceR[i] = forceX0 + (forceX1 - forceX0)*i/(nR_force-1);
      }
      for(int i=0;i<nZ_force;i++)
      {
        forceZ[i] = forceZ0 + (forceZ1 - forceZ0)*i/(nZ_force-1);
      }
      float Btotal = 0.0;
      for(int i=0;i<nR_force;i++)
      {
        for(int j=0;j<nZ_force;j++)
        {
            positionHistoryX[i][tt/subSampleFac] = particleArray->xprevious[i];
            positionHistoryY[i][tt/subSampleFac] = particleArray->yprevious[i];
            positionHistoryZ[i][tt/subSampleFac] = particleArray->zprevious[i];
            velocityHistoryX[i][tt/subSampleFac] = particleArray->vx[i];
            velocityHistoryY[i][tt/subSampleFac] = particleArray->vy[i];
            velocityHistoryZ[i][tt/subSampleFac] = particleArray->vz[i];
            chargeHistory[i][tt/subSampleFac] = particleArray->charge[i];
        }
}
#endif
#endif
//#else
/*
#if USE_BOOST
cpu_times moveTime0 = timer.elapsed();
#endif
        std::for_each(particleArray.begin(), particleArray.end(),
                move_boris(dt,boundaries.data(),nLines, 
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front(),
                    nR_PreSheathEfield,nZ_PreSheathEfield, 
                    &preSheathEGridr.front(),&preSheathEGridz.front(),
                    &PSEr.front(),&PSEz.front(),&PSEt.front()));
#if USE_BOOST
cpu_times moveTime1 = timer.elapsed();
moveTime = moveTime + (moveTime1.wall - moveTime0.wall);
cpu_times geomTime0 = timer.elapsed();
#endif
    std::for_each(particleArray.begin(), particleArray.end(),
            geometry_check(nLines,boundaries.data(),dt,tt) );
#if USE_BOOST
cpu_times geomTime1 = timer.elapsed();
geomCheckTime = geomCheckTime + (geomTime1.wall - geomTime0.wall);
#endif
#if USEIONIZATION > 0
#if USE_BOOST
cpu_times ionizTime0 = timer.elapsed();
#endif
std::cout << "Flow vNs "<< testFlowVec[0] << " " <<testFlowVec[1] << " " << testFlowVec[2]  << std::endl;
    std::cout << "Starting main loop" << particleArray->xprevious[0] << std::endl;
    std::cout << "Starting main loop" << particleArray->xprevious[1] << std::endl;
//Main time loop
    #if __CUDACC__
      cudaDeviceSynchronize();
    #endif
     //int nDevices=0;
     //nDevices = omp_get_num_threads();
     //    unsigned int cpu_thread_id = omp_get_thread_num();
     //    unsigned int num_cpu_threads = omp_get_num_threads();
     //    printf("Number of CPU threads %d (ID %d)\n", cpu_thread_id, num_cpu_threads);
    #if USE_OPENMP
    //for(int device=0;device <1;device++)
//{ cudaSetDevice(device); 
    //int nDevices=32;
    //std::cout << "nDevices " << nDevices << std::endl;
     //omp_set_num_threads(nDevices);  // create as many CPU threads as there are CUDA devices
    //std::cout << "nDevices " << nDevices << std::endl;
     int nDevices=0;
     #pragma omp parallel
     {
        nDevices = omp_get_num_threads();
        //tid = omp_get_thread_num();
         unsigned int cpu_thread_id = omp_get_thread_num();
         unsigned int num_cpu_threads = omp_get_num_threads();
         //int gpu_id = -1;
             //cudaSetDevice(cpu_thread_id);
        //cudaSetDevice(cpu_thread_id % nDevices);        // "% num_gpus" allows more CPU threads than GPU devices
        //if(cpu_thread_id ==0 ) cudaSetDevice(0); 
        //if(cpu_thread_id ==1 ) cudaSetDevice(3); 
        //cudaGetDevice(&gpu_id);
         printf("CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads);
#if USE_MPI > 0 
    printf("Hello world from processor %s, rank %d"
           " out of %d processors and cpu_thread_id %i \n",
                      processor_name, world_rank, world_size,cpu_thread_id);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
#if USERECOMBINATION > 0
    std::for_each(particleArray.begin(), particleArray.end(), 
            recombine(dt) );
#endif
sim::Array<int> tmpInt(1,1),tmpInt2(1,1);
   //int nN=10000;
   //thrust::host_vector<int> h_vec(nN); 
   //thrust::generate(h_vec.begin(), h_vec.end(), rand);
   //// transfer data to the device
   //thrust::device_vector<int> d_vec = h_vec;
   //float *d_vec2;
   //cudaMallocManaged(&d_vec2, 1000*sizeof(float));
   //std::cout << "created d_vec and cmalloc, starting init " << std::endl; 
   //for(int k=0;k<1000;k++)
   //{   //std::cout << "k " << k << std::endl;
   //    d_vec2[k] = 1.0f;
   //}
   //for(int k=0;k<1000;k++)
   //{   //std::cout << "k " << k << std::endl;
   //    //d_vec2[k] = 1.0f;
   //thrust::sort(thrust::device,d_vec.begin()+world_rank*nN/world_size, d_vec.begin()+ (world_rank+1)*nN/world_size-1); // sort data on the device 
   //}
   //// transfer data back to host
   //thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
    for(tt; tt< nT; tt++)
    {
     dev_tt[0] = tt;
     std::cout << "beginning of loop tt " << tt << std::endl;
      #if USE_SORT > 0
        thrust::for_each(thrust::device,tmpInt.begin(),tmpInt.end(),sort0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif

      #if PARTICLE_TRACKS >0
     std::cout << "particle tracks" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],history0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif
      //std::cout << " world rank pstart nactive " << world_rank << " " << pStartIndx[world_rank] << "  " << nActiveParticlesOnRank[world_rank] << std::endl;
      //thrust::for_each(thrust::device,particleBegin,particleOne,
      //     test_routinePp(particleArray));
     std::cout << "boris" << std::endl;

      thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],
                move_boris0);
      #ifdef __CUDACC__
        cudaDeviceSynchronize();
      #endif
     std::cout << "check geom" << std::endl;
      thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],
                    geometry_check0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif

      #if SPECTROSCOPY > 0
     std::cout << "spec" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],spec_bin0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif  

      #if USEIONIZATION > 0
     std::cout << "ioni" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],
                ionize0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif

      #if USERECOMBINATION > 0
     std::cout << "rec" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],recombine0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif

      #if USEPERPDIFFUSION > 0
     std::cout << "diff" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],crossFieldDiffusion0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
            
     std::cout << "diff geom check" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],geometry_check0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif

      #if USECOULOMBCOLLISIONS > 0
     std::cout << "coll" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],coulombCollisions0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif
      
      #if USETHERMALFORCE > 0
     std::cout << "therm" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],thermalForce0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif

      #if USESURFACEMODEL > 0
     std::cout << "surf" << std::endl;
        thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],reflection0);
      #ifdef __CUDACC__
        cudaThreadSynchronize();
      #endif
      #endif        

    }
   #if PARTICLE_TRACKS >0
     tt = nT;
     std::cout << " tt for final history " << tt << std::endl;
     thrust::for_each(thrust::device,particleBegin+pStartIndx[world_rank],particleBegin+pStartIndx[world_rank]+nActiveParticlesOnRank[world_rank],
      history0);
   #endif

#if USE_OPENMP
}
#endif
*/
//#endif
    }
// Ensure that all time step loop GPU kernels are complete before proceeding
    #ifdef __CUDACC__
        cudaDeviceSynchronize();
    #endif

    auto finish_clock = Time::now();
    fsec fs = finish_clock - start_clock;
    printf("Time taken          is %6.3f (secs) \n", fs.count());
    printf("Time taken per step is %6.3f (secs) \n", fs.count() / (float) nT);
#if USE_BOOST
    cpu_times ionizeTimeGPU = timer.elapsed();
    std::cout << "Particle Moving Time: " << ionizeTimeGPU.wall*1e-9 << '\n';
#endif
    /*
for(int i=0; i<nP ; i++)
{
    std::cout << "particle " << i << " first rnd# " << 
        particleArray->test[i] << " and x " << particleArray->xprevious[i] << 
         " hitwall " << particleArray->hitWall[i] << 
         " trans " << particleArray->transitTime[i] << std::endl;
}
*/
    std::cout << "transit time counting "<< nP << " " << particleArray->x[0] <<  std::endl;
    //float tmp202 =0.0;
#if USE_CUDA
    cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
    sim::Array<float> xGather(nP,0.0);
    sim::Array<float> test0Gather(nP,0.0);
    sim::Array<float> test1Gather(nP,0.0);
    sim::Array<float> yGather(nP,0.0);
    sim::Array<float> zGather(nP,0.0);
    sim::Array<float> vGather(nP,0.0);
    sim::Array<float> vxGather(nP,0.0);
    sim::Array<float> vyGather(nP,0.0);
    sim::Array<float> vzGather(nP,0.0);
    sim::Array<float> hitWallGather(nP,0.0);
    sim::Array<float> weightGather(nP,0.0);
    sim::Array<float> chargeGather(nP,0.0);
    sim::Array<float> firstIonizationTGather(nP,0.0);
    sim::Array<float> firstIonizationZGather(nP,0.0);
    //float *x_gather = NULL;
    //if (world_rank == 0) {
    //      x_gather = malloc(sizeof(float)*nP);
    //}
    std::cout << "reached gather barrier" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "started gather" << std::endl;
    MPI_Gather(&particleArray->x[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &xGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->y[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &yGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->z[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &zGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->v[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &vGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->vx[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &vxGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->vy[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &vyGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->vz[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &vzGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->hitWall[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &hitWallGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->weight[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &weightGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->charge[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &chargeGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->firstIonizationT[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &firstIonizationTGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->firstIonizationZ[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &firstIonizationZGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->test0[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &test0Gather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&particleArray->test1[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &test1Gather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);
    std::cout << "wating after gathers" << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "passed barrier after gather" << std::endl;
#if PARTICLE_TRACKS >0
   
    std::vector<float> exampleArray(4,0.0);
    std::vector<float> exampleArrayGather(4,0.0);
    if(world_rank ==0)
    {
      exampleArray[0]=1;
      exampleArray[1]=1;
    }
    if(world_rank ==1)
    {
      exampleArray[2]=2;
      exampleArray[3]=2;
    }
    std::vector<int> exCount(2,2),exDispl(2,0);
    exDispl[0] = 0;
    exDispl[1] = 2;
    const int* exdispl=&exDispl[0];
    const int* excount = &exCount[0];

    //MPI_Gatherv(&exampleArray[exDispl[world_rank]],2,MPI_FLOAT,&exampleArrayGather[0],excount,exdispl,MPI_FLOAT,0,MPI_COMM_WORLD);

    //for(int i=0;i<4;i++)
    //{
    //  std::cout << "rank " << world_rank << " val " << exampleArrayGather[i] << std::endl; 
    //}

MPI_Barrier(MPI_COMM_WORLD);

    //for(int i=pDisplacement[world_rank];i<pDisplacement[world_rank]+pHistPerNode[world_rank];i++)
    //{
    //  std::cout << "Rank i "<< i << " "  << world_rank << "z " << positionHistoryZ[i] << std::endl;
    //}
    //std::cout << "starting particle tracks gather "<< world_rank<< " pstart "<< pStartIndx[world_rank] << "nhist " << nHistoriesPerParticle << std::endl;
    //std::cout << "start gather 2 "<< world_rank<< " nppr "<< nPPerRank[world_rank] << "nhist " << nHistoriesPerParticle << std::endl;
    MPI_Gatherv(&positionHistoryX[pDisplacement[world_rank]], pHistPerNode[world_rank],MPI_FLOAT,&positionHistoryXgather[0],phpn,displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&positionHistoryY[pDisplacement[world_rank]], pHistPerNode[world_rank],MPI_FLOAT,&positionHistoryYgather[0],phpn,displ, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&positionHistoryZ[pDisplacement[world_rank]], pHistPerNode[world_rank], MPI_FLOAT, &positionHistoryZgather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&velocityHistory[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &velocityHistorygather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&velocityHistoryX[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &velocityHistoryXgather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&velocityHistoryY[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &velocityHistoryYgather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&velocityHistoryZ[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &velocityHistoryZgather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&chargeHistory[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &chargeHistoryGather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&weightHistory[pStartIndx[world_rank]*nHistoriesPerParticle], nPPerRank[world_rank]*nHistoriesPerParticle, MPI_FLOAT, &weightHistoryGather[0], phpn,displ,MPI_FLOAT, 0, MPI_COMM_WORLD);
    std::cout << "at barrier tracks gather" << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "finished particle tracks gather" << std::endl;
    //if(world_rank ==0)
    //{
    //for(int i=0;i<401;i++)
    //{
    //  std::cout << "Rank " << world_rank << "z " << positionHistoryZgather[i] << std::endl;
    //}
    //}
#endif

#if SPECTROSCOPY > 0
MPI_Barrier(MPI_COMM_WORLD);
std::cout <<"Starting spectroscopy reduce " << world_rank<< std::endl;
MPI_Reduce(&net_Bins[0], &net_BinsTotal[0], (nBins+1)*net_nX*net_nY*net_nZ, 
           MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
std::cout <<"Done with spectroscopy reduce " << world_rank<< std::endl;
MPI_Barrier(MPI_COMM_WORLD);
#endif
#if (USESURFACEMODEL > 0 || FLUX_EA > 0)
//MPI_Barrier(MPI_COMM_WORLD);
std::cout <<"Starting surface reduce " << world_rank<< std::endl;
//for(int i=0;i<nSurfaces;i++) std::cout << surfaces->grossDeposition[i]<<std::endl;
MPI_Reduce(&surfaces->grossDeposition[0], &grossDeposition[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->grossErosion[0], &grossErosion[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sumWeightStrike[0], &sumWeightStrike[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->aveSputtYld[0], &aveSputtYld[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sputtYldCount[0], &sputtYldCount[0],nSurfaces, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sumParticlesStrike[0], &sumParticlesStrike[0],nSurfaces, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Reduce(&surfaces->energyDistribution[0], &energyDistribution[0],nSurfaces*nEdist*nAdist, 
        MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sputtDistribution[0], &sputtDistribution[0],nSurfaces*nEdist*nAdist, 
        MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Reduce(&surfaces->reflDistribution[0], &reflDistribution[0],nSurfaces*nEdist*nAdist, 
        MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    if(world_rank == 0)
{
    auto MPIfinish_clock = Time::now();
    fsec fsmpi = MPIfinish_clock - finish_clock;
    printf("Time taken for mpi reduction          is %6.3f (secs) \n", fsmpi.count());
}
//    tmp202 =  particleArray->vx[0];
    //std::cout << "memory access hitwall " 
    //<< particleArray->xprevious[0] << std::endl;
    //std::cout << "transit time counting " << std::endl;
  #if USE_MPI > 0
    if(world_rank == 0)
      {
  #endif
        int totalHitWall=0;
        for(int i=0;i<nP;i++)
        {
          if(particleArray->hitWall[i] > 0.0) totalHitWall++;
        }
        std::cout << "Number and percent of particles that hit wall " << 
        totalHitWall << " " << totalHitWall*1.0/(nP*1.0) << std::endl;
        #if USE3DTETGEOM > 0
          float meanTransitTime0 = 0.0;
    /*
    for (int i=0; i<nP; i++)
    {
        if(particleArray->hitWall[i] == 1.0)
        {
            meanTransitTime0 = meanTransitTime0 + particleArray->transitTime[i];
        }
    }
meanTransitTime0 = meanTransitTime0/nP;
std::cout << " mean transit time " << meanTransitTime0 << std::endl;
    int max_boundary = 0;
    float max_impacts = 0.0;
    int max_boundary1 = 0;
    float max_impacts1 = 0.0;
    float* impacts = new float[nLines];
    for (int i=0; i<nLines; i++)
    {
        impacts[i] = boundaries[i].impacts;
        if (boundaries[i].impacts > max_impacts)
        {
            max_impacts = boundaries[i].impacts;
            max_boundary = i;
        }
    }


std::cout << "maximum boundary " << max_boundary << std::endl;
std::cout << "number of counts " << max_impacts << std::endl;
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
    float* impacts = new float[nLines];
    for (int i=0; i<nLines; i++)
    {
        impacts[i] = boundaries[i].impacts;
    }
#endif
#if PARTICLE_SOURCE == 1
int ring1 = 0;
int ring2 = 0;
int noWall = 0;
float meanTransitTime = 0.0;

for(int i=0; i<nP ; i++)
{
	if(particleArray->wallIndex[i] == boundaryIndex_ImpurityLaunch[0])
	{
		ring1++;
	}
	else if(particleArray->wallIndex[i] == boundaryIndex_ImpurityLaunch[1])
	{
		ring2++;
	}
	
	if(particleArray->wallIndex[i] == 0)
	{
		noWall++;
	}
	
	meanTransitTime = meanTransitTime + particleArray->transitTime[i];
	
} 
meanTransitTime = meanTransitTime/(nP-noWall);
std::cout << "Number of impurity particles deposited on ring 1 " << ring1 << std::endl;
std::cout << "Number of impurity particles deposited on ring 2 " << ring2 << std::endl;
std::cout << "Number of impurity particles not deposited " << noWall << std::endl;
std::cout << "Mean transit time of deposited particles " << meanTransitTime << std::endl;
#endif
    ofstream outfile2;
    outfile2.open ("positions.m");
    for(int i=1 ; i<=nP ; i++)
      {
        outfile2 << "Pos( " << i<< ",:) = [ " ;
        outfile2 << particleArray->x[i-1] << " " << particleArray->y[i-1] 
            << " " << particleArray->z[i-1] << " ];" << std::endl;
      }
       outfile2.close();
// Write netCDF output for positions
for (int i=0; i<nP; i++)
{
    finalPosX[i] = particleArray->xprevious[i];
    finalPosY[i] = particleArray->yprevious[i];
    finalPosZ[i] = particleArray->zprevious[i];
    finalVx[i] =   particleArray->vx[i];
    finalVy[i] =   particleArray->vy[i];
    finalVz[i] =   particleArray->vz[i];
    transitTime[i] = particleArray->transitTime[i];
}
NcFile ncFile0("positions.nc", NcFile::replace);
NcDim nc_nP0 = ncFile0.addDim("nP",nP);
vector<NcDim> dims0;
dims0.push_back(nc_nP0);

NcVar nc_x0 = ncFile0.addVar("x",ncDouble,dims0);
NcVar nc_y0 = ncFile0.addVar("y",ncDouble,dims0);
NcVar nc_z0 = ncFile0.addVar("z",ncDouble,dims0);
NcVar nc_vx0 = ncFile0.addVar("vx",ncDouble,dims0);
NcVar nc_vy0 = ncFile0.addVar("vy",ncDouble,dims0);
NcVar nc_vz0 = ncFile0.addVar("vz",ncDouble,dims0);
NcVar nc_trans0 = ncFile0.addVar("transitTime",ncDouble,dims0);

nc_x0.putVar(finalPosX);
nc_y0.putVar(finalPosY);
nc_z0.putVar(finalPosZ);
nc_vx0.putVar(finalVx);
nc_vy0.putVar(finalVy);
nc_vz0.putVar(finalVz);
nc_trans0.putVar(transitTime);

NcFile ncFile1("surface.nc", NcFile::replace);
NcDim nc_nLines = ncFile1.addDim("nLines",nLines);
vector<NcDim> dims1;
dims1.push_back(nc_nLines);

vector<NcDim> dimsSurfE;
dimsSurfE.push_back(nc_nLines);
NcDim nc_nEnergies = ncFile1.addDim("nEnergies",nEdist);
NcDim nc_nAngles = ncFile1.addDim("nAngles",nAdist);
dimsSurfE.push_back(nc_nAngles);
dimsSurfE.push_back(nc_nEnergies);
NcVar nc_grossDep = ncFile1.addVar("grossDeposition",ncFloat,nc_nLines);
NcVar nc_grossEro = ncFile1.addVar("grossErosion",ncFloat,nc_nLines);
NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl",ncFloat,nc_nLines);
NcVar nc_spylCounts = ncFile1.addVar("spylCounts",ncInt,nc_nLines);
NcVar nc_surfNum = ncFile1.addVar("surfaceNumber",ncInt,nc_nLines);
NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike",ncInt,nc_nLines);
NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike",ncFloat,nc_nLines);
nc_grossDep.putVar(&grossDeposition[0]);
nc_surfNum.putVar(&surfaceNumbers[0]);
nc_grossEro.putVar(&grossErosion[0]);
nc_aveSpyl.putVar(&aveSputtYld[0]);
nc_spylCounts.putVar(&sputtYldCount[0]);
nc_sumParticlesStrike.putVar(&sumParticlesStrike[0]);
nc_sumWeightStrike.putVar(&sumWeightStrike[0]);
//NcVar nc_surfImpacts = ncFile1.addVar("impacts",ncFloat,dims1);
//NcVar nc_surfRedeposit = ncFile1.addVar("redeposit",ncFloat,dims1);
//NcVar nc_surfStartingParticles = ncFile1.addVar("startingParticles",ncFloat,dims1);
//NcVar nc_surfZ = ncFile1.addVar("Z",ncFloat,dims1);
NcVar nc_surfEDist = ncFile1.addVar("surfEDist",ncFloat,dimsSurfE);
NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist",ncFloat,dimsSurfE);
NcVar nc_surfSputtDist = ncFile1.addVar("surfSputtDist",ncFloat,dimsSurfE);
//nc_surfImpacts.putVar(impacts);
//#if USE3DTETGEOM > 0
//nc_surfRedeposit.putVar(redeposit);
//#endif
//nc_surfStartingParticles.putVar(startingParticles);
//nc_surfZ.putVar(surfZ);
std::cout << "writing energy distribution file " << std::endl;
nc_surfEDist.putVar(&energyDistribution[0]);
nc_surfReflDist.putVar(&reflDistribution[0]);
nc_surfSputtDist.putVar(&sputtDistribution[0]);
//NcVar nc_surfEDistGrid = ncFile1.addVar("gridE",ncDouble,nc_nEnergies);
//nc_surfEDistGrid.putVar(&surfaces->gridE[0]);
//NcVar nc_surfADistGrid = ncFile1.addVar("gridA",ncDouble,nc_nAngles);
//nc_surfADistGrid.putVar(&surfaces->gridA[0]);
ncFile1.close();
#endif
#if PARTICLE_TRACKS > 0

// Write netCDF output for histories
NcFile ncFile_hist("history.nc", NcFile::replace);
NcDim nc_nT = ncFile_hist.addDim("nT",nT/subSampleFac);
NcDim nc_nP = ncFile_hist.addDim("nP",nP);
vector<NcDim> dims_hist;
#if USE_CUDA
NcDim nc_nPnT = ncFile_hist.addDim("nPnT",nP*nT/subSampleFac);
dims_hist.push_back(nc_nPnT);
#else
dims_hist.push_back(nc_nP);
dims_hist.push_back(nc_nT);
#endif
NcVar nc_x = ncFile_hist.addVar("x",ncDouble,dims_hist);
NcVar nc_y = ncFile_hist.addVar("y",ncDouble,dims_hist);
NcVar nc_z = ncFile_hist.addVar("z",ncDouble,dims_hist);

NcVar nc_v = ncFile_hist.addVar("v",ncDouble,dims_hist);
NcVar nc_vx = ncFile_hist.addVar("vx",ncDouble,dims_hist);
NcVar nc_vy = ncFile_hist.addVar("vy",ncDouble,dims_hist);
NcVar nc_vz = ncFile_hist.addVar("vz",ncDouble,dims_hist);

NcVar nc_charge = ncFile_hist.addVar("charge",ncDouble,dims_hist);
NcVar nc_weight = ncFile_hist.addVar("weight",ncDouble,dims_hist);
#if USE_MPI > 0
    //if(world_rank ==0)
    //{
    //for(int i=0;i<401;i++)
    //{
    //  std::cout << "Rank " << world_rank << "z " << positionHistoryZgather[i] << std::endl;
    //}
    //}
std::cout << "printing gathered data" << std::endl;
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
nc_x.putVar(positionHistoryX[0]);
nc_y.putVar(positionHistoryY[0]);
nc_z.putVar(positionHistoryZ[0]);

nc_vx.putVar(velocityHistoryX[0]);
nc_vy.putVar(velocityHistoryY[0]);
nc_vz.putVar(velocityHistoryZ[0]);

nc_charge.putVar(chargeHistory[0]);
#endif
#endif
#if SPECTROSCOPY > 0
// Write netCDF output for density data
NcFile ncFile("spec.nc", NcFile::replace);
NcDim nc_nBins = ncFile.addDim("nBins",nBins+1);
NcDim nc_nR = ncFile.addDim("nR",net_nX);
NcDim nc_nZ = ncFile.addDim("nZ",net_nZ);

vector<NcDim> dims;
dims.push_back(nc_nBins);
dims.push_back(nc_nZ);
dims.push_back(nc_nR);

NcVar nc_n = ncFile.addVar("n",ncDouble,dims);
float *binPointer = &net_Bins[0];
nc_n.putVar(binPointer);
#endif
#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif
#if USE_BOOST
    cpu_times copyToHostTime = timer.elapsed();

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "Copy to host, bin and output time: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';
    std::cout << "Total ODE integration time: " << moveTime*1e-9 << '\n';
    std::cout << "Total geometry checking time: " << geomCheckTime*1e-9 << '\n';
    std::cout << "Total ionization time: " << ionizTime*1e-9 << '\n';
#endif
#if USE_MPI > 0
    std::cout << "sqrt 0 " << sqrt(-0.0) << std::endl;
    for(int i=0;i<100;i++)
{
    //std::cout << "particle hitwall and Ez " << particleArray->hitWall[i] << " " << particleArray->test[i] << " "<< test0Gather[i] << " " << test1Gather[i]<< " "<<
    //    particleArray->test2[i] << " " << particleArray->test3[i] << " " << particleArray->test4[i] << 
    //    " " << particleArray->distTraveled[i] << std::endl;
}
    for(int i=0;i<100;i++)
{
    //std::cout << "particle ionization z and t " << firstIonizationZGather[i] << " " << firstIonizationTGather[i] << " " << xGather[i] << " " << 
      //vxGather[i] << " " << chargeGather[i] << std::endl;
}
#endif
//    for(int i=0;i<100;i++)
//{
//    std::cout << "reflected/sputtered energy " << particleArray->newVelocity[i]   << std::endl;
//}
//#endif
#if USE_MPI > 0
    }
#endif

#ifdef __CUDACC__
    cudaError_t err = cudaDeviceReset();
//cudaProfilerStop();
#endif
#if USE_MPI > 0
    // Finalize the MPI environment.
    MPI_Finalize();
#endif
  if(world_rank == 0)
  {
    auto GITRfinish_clock = Time::now();
    fsec fstotal = GITRfinish_clock - GITRstart_clock;
    printf("Total runtime for GITR is %6.3f (secs) \n", fstotal.count());
  }
//#endif
    return 0;
    }
