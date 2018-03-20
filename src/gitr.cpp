#include <iostream>
#include <stdio.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <vector>
#include <netcdf>
#include "ncFile.h"
#include "utils.h"
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
#include "Particles.h"
#include "Boundary.h"
#include "Surfaces.h"
#include "curandInitialize.h"
#include "spectroscopy.h"
#include "fieldLineTrace.h"
#include "history.h"
#include "io.hpp"
#include "hashGeom.h"
#include "hashGeomSheath.h"
#include "testRoutine.h"
#include "testRoutineCuda.h"
#include "boundaryInit.h"
#include "array.h"
#include "ompPrint.h"

#if USE_BOOST
    #include <boost/timer/timer.hpp>
    #include "boost/filesystem.hpp"
    #include <boost/random.hpp>
#endif

#ifdef __CUDACC__
    #include <curand.h>
    #include <curand_kernel.h>
#endif

#if USE_MPI 
#include <mpi.h>
#endif

#if USE_OPENMP
#include <omp.h>
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

int main(int argc, char **argv)
{
 #if USE_MPI > 0
        int ppn = 4;
        int np = 1;
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);
    
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
    // Finalize the MPI environment.
    //MPI_Finalize();
  #else
    int world_rank=0;
    int world_size=1;
  #endif
  //Prepare config files for import
  Config cfg,cfg_geom;
  std::string input_path = "input/";
  //Parse and read input file
  std::cout << "Open configuration file input/gitrInput.cfg " << std::endl;
  importLibConfig(cfg,input_path+"gitrInput.cfg");
  
  // Parse and read geometry file
  std::string geomFile; 
  getVariable(cfg,"geometry.fileString",geomFile);
  std::cout << "Open geometry file " << input_path+geomFile << std::endl; 
  importLibConfig(cfg_geom,input_path+geomFile);

  std::cout << "Successfully staged input and geometry file " << std::endl;
  
  //check binary compatibility with input file
  #if CHECK_COMPATIBILITY>0
    checkFlags(cfg); 
  #endif
  
  // show memory usage of GPU
  #if USE_CUDA 
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t    cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  
    if(cudaSuccess != cuda_status )
    {
  
       printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
       exit(1);
    }
  
    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;
    
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
      used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0); 
    int nDevices;
    int nThreads;
    cudaGetDeviceCount(&nDevices);
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
  #endif 
  #if USE_BOOST
    //Output
    boost::filesystem::path dir("output"); 
    if(!(boost::filesystem::exists(dir)))
    {
      // std::cout<<"Doesn't Exists"<<std::endl;
      if (boost::filesystem::create_directory(dir))
      {
         //std::cout << " Successfully Created " << std::endl;
      }
    }
  #endif
  // Background species info
  float background_Z,background_amu;
  getVariable(cfg,"backgroundPlasmaProfiles.Z",background_Z);
  getVariable(cfg,"backgroundPlasmaProfiles.amu",background_amu);

  //Bfield initialization
  int nR_Bfield = 1;
  int nY_Bfield = 1;
  int nZ_Bfield = 1;
  int n_Bfield = 1;
  std::string bfieldCfg = "backgroundPlasmaProfiles.Bfield.";
  #if BFIELD_INTERP > 0
    std::string bfieldFile;
    getVariable(cfg,bfieldCfg+"fileString",bfieldFile);
    nR_Bfield = getDimFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridNrString");
  #endif
  #if BFIELD_INTERP > 1
    nZ_Bfield = getDimFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridNzString");
  #endif
  #if BFIELD_INTERP > 2
    nY_Bfield = getDimFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridNyString");
  #endif
  sim::Array<float> bfieldGridr(nR_Bfield),bfieldGridy(nY_Bfield),bfieldGridz(nZ_Bfield);
  n_Bfield = nR_Bfield*nY_Bfield*nZ_Bfield;
  sim::Array<float> br(n_Bfield),by(n_Bfield),bz(n_Bfield);
  #if BFIELD_INTERP == 0
    getVariable(cfg,bfieldCfg+"r",br[0]);
    getVariable(cfg,bfieldCfg+"y",by[0]);
    getVariable(cfg,bfieldCfg+"z",bz[0]);
  #else
    getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridRString",bfieldGridr);
    #if BFIELD_INTERP > 1
      getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridZString",bfieldGridz);
    #endif
    #if BFIELD_INTERP > 2
      getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"gridYString",bfieldGridy);
    #endif

    getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"rString",br);
    getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"yString",by);
    getVarFromFile(cfg,input_path+bfieldFile,bfieldCfg,"zString",bz);
  #endif 
  float Btest[3] = {0.0f}; 
  interp2dVector(&Btest[0],0.0,0.0,0.0,nR_Bfield,
                    nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
  std::cout << "Bfield at 000 "<< Btest[0] << " " << Btest[1] << " " << Btest[2] << std::endl; 
  interp2dVector(&Btest[0],0.0,0.0,0.2,nR_Bfield,
                    nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
  std::cout << "Bfield at 000.2 "<< Btest[0] << " " << Btest[1] << " " << Btest[2] << std::endl; 
  std::cout << "Finished Bfield import" << std::endl;
#if USE_MPI > 0 
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Bcast(br.data(), nR_Bfield,MPI_FLOAT,0,MPI_COMM_WORLD);
  //MPI_Bcast(bz.data(), nZ_Bfield,MPI_FLOAT,0,MPI_COMM_WORLD);
  //printf("Hello world from processor %s, rank %d"
  //                   " out of %d processors bfield value %f \n",
  //                                         processor_name, world_rank, world_size,br[0]);
#endif
  std::string profiles_folder = "output/profiles";  
  
  //Geometry Definition
  Setting& geom = cfg_geom.lookup("geom");
  int nLines = geom["x1"].getLength();
  int nSurfaces=0;
  //int nMaterials = geom["nMaterials"];
  std::cout << "Number of Geometric Objects To Load: " << nLines << std::endl;
  
  sim::Array<Boundary> boundaries(nLines+1,Boundary());
  nSurfaces = importGeometry(cfg_geom, boundaries);
  std::cout << "Starting Boundary Init..." << std::endl;
  std::cout << " y1 and y2 " << boundaries[nLines].y1 << " " << boundaries[nLines].y2 << std::endl;
  float biasPotential = 0.0;
  #if BIASED_SURFACE > 0
    getVariable(cfg,"backgroundPlasmaProfiles.biasPotential",biasPotential);
  #endif
  //create Surface data structures
  int nEdist = 1;
  float E0dist = 0.0;
  float Edist = 0.0;
  int nAdist = 1;
  float A0dist = 0.0;
  float Adist = 0.0;
  #if FLUX_EA > 0
    getVariable(cfg,"surfaces.flux.nE",nEdist);
    getVariable(cfg,"surfaces.flux.E0",E0dist);
    getVariable(cfg,"surfaces.flux.E",Edist);

    getVariable(cfg,"surfaces.flux.nA",nAdist);
    getVariable(cfg,"surfaces.flux.A0",A0dist);
    getVariable(cfg,"surfaces.flux.A",Adist);
    std::cout << "dist nE E0 E nA A0 A" << nEdist << " " << E0dist << " " 
       << Edist << " " << nAdist << " " << A0dist << " " << Adist << std::endl;
    std::cout << "nLines before surfaces " << nLines << std::endl;
  #endif
  auto surfaces = new Surfaces(nSurfaces,nEdist,nAdist);
  std::cout << "nLines after surfaces " << sizeof(surfaces) << std::endl;
  //surfaces->setSurface(nEdist, E0dist,Edist,nAdist ,A0dist ,Adist);
  std::cout << "surfaces sizeof" << sizeof(surfaces) << std::endl;

  std::cout << "surface stuff " << surfaces->nE << " " << surfaces->E0 << " " << surfaces->E << " " << surfaces->dE <<  std::endl;
  std::cout << "surface stuff " << surfaces->nA << " " << surfaces->A0 << " " << surfaces->A << " " << surfaces->dA <<  std::endl;
  sim::Array<float> grossDeposition(nSurfaces,0.0);
  sim::Array<float> grossErosion(nSurfaces,0.0);
  sim::Array<float> sumWeightStrike(nSurfaces,0.0);
  sim::Array<float> energyDistribution(nSurfaces*nEdist*nAdist,0.0);
  sim::Array<float> aveSputtYld(nSurfaces,0.0);
  sim::Array<int> sputtYldCount(nSurfaces,0);
  sim::Array<int> sumParticlesStrike(nSurfaces,0);

  int nR_closeGeom = 1;
  int nY_closeGeom = 1;
  int nZ_closeGeom = 1;
  int n_closeGeomElements = 1;
  int nGeomHash = 1;
  std::string geomHashCfg = "geometry_hash.";
  #if GEOM_HASH == 1
    getVariable(cfg,geomHashCfg+"nR_closeGeom",nR_closeGeom);
    getVariable(cfg,geomHashCfg+"nZ_closeGeom",nZ_closeGeom);
    getVariable(cfg,geomHashCfg+"n_closeGeomElements",n_closeGeomElements);
    nGeomHash = nR_closeGeom*nZ_closeGeom*n_closeGeomElements;
    #if USE3DTETGEOM > 0
      getVariable(cfg,geomHashCfg+"nY_closeGeom",nY_closeGeom);
      nGeomHash = nY_closeGeom*nGeomHash;
    #endif
  #endif

  #if GEOM_HASH > 1
    std::string hashFile;
    getVariable(cfg,geomHashCfg+"fileString",hashFile);
    nR_closeGeom = getDimFromFile(cfg,input_path+hashFile,geomHashCfg,"gridNrString");
    nZ_closeGeom = getDimFromFile(cfg,input_path+hashFile,geomHashCfg,"gridNzString");
    n_closeGeomElements = getDimFromFile(cfg,input_path+hashFile,geomHashCfg,"nearestNelementsString");
    nGeomHash = nR_closeGeom*nZ_closeGeom*n_closeGeomElements;
    #if USE3DTETGEOM > 0
      nY_closeGeom = getDimFromFile(cfg,input_path+hashFile,geomHashCfg,"gridNyString");
      nGeomHash = nY_closeGeom*nGeomHash;
    #else
    #endif
  #endif
  sim::Array<float> closeGeomGridr(nR_closeGeom), closeGeomGridy(nY_closeGeom),
      closeGeomGridz(nZ_closeGeom);
  sim::Array<int> closeGeom(nGeomHash,0);
  #if GEOM_HASH == 1
    float hashX0,hashX1,hashY0,hashY1,hashZ0,hashZ1;
    getVariable(cfg,geomHashCfg+"hashX0",hashX0);
    getVariable(cfg,geomHashCfg+"hashX1",hashX1);
    getVariable(cfg,geomHashCfg+"hashZ0",hashZ0);
    getVariable(cfg,geomHashCfg+"hashZ1",hashZ1);
    #if USE3DTETGEOM > 0
      getVariable(cfg,geomHashCfg+"hashY0",hashY0);
      getVariable(cfg,geomHashCfg+"hashY1",hashY1);
    #endif
    
    for(int i=0; i<nR_closeGeom; i++)
    {  closeGeomGridr[i] = (hashX1 - hashX0)*i/(nR_closeGeom - 1)+ hashX0;}
    for(int j=0; j<nY_closeGeom; j++)
    {  closeGeomGridy[j] = (hashY1 - hashY0)*j/(nY_closeGeom - 1)+ hashY0;}
    for(int k=0; k<nZ_closeGeom; k++)
    {  closeGeomGridz[k] = (hashZ1 - hashZ0)*k/(nZ_closeGeom - 1)+ hashZ0;}
   std::cout << "about to create iterator1 " << std::endl; 
    thrust::counting_iterator<std::size_t> lines0(0);  
   std::cout << "iterator2 " << std::endl; 
    thrust::counting_iterator<std::size_t> lines1(nR_closeGeom*nY_closeGeom*nZ_closeGeom);
    int nHashMeshPointsPerProcess=ceil(nR_closeGeom*nY_closeGeom*nZ_closeGeom/world_size);
   std::cout << "nHashMeshPointsPerProcess "<< nHashMeshPointsPerProcess << std::endl; 
   std::vector<int> hashMeshIncrements(world_size);
   for(int j=0;j<world_size-1;j++)
   {
       hashMeshIncrements[j] = nHashMeshPointsPerProcess;
   }
   hashMeshIncrements[world_size-1]=nR_closeGeom*nY_closeGeom*nZ_closeGeom - (world_size-1)*nHashMeshPointsPerProcess;
   std::cout << "minDist1 "<< nGeomHash << std::endl; 
    sim::Array<float> minDist1(n_closeGeomElements,1e6);
    std::cout << "Generating geometry hash" << sizeof(int) << " bytes per int, " 
        << nGeomHash << " for the entire hash " <<  std::endl;

#if USE_CUDA >0
     cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  
    if(cudaSuccess != cuda_status )
    {
  
       printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
       exit(1);
    }
  
    free_db = (double)free_byte ;
    total_db = (double)total_byte ;
    used_db = total_db - free_db ;
    
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
      used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0); 
  #endif
    //for(int i=0; i<nZ_closeGeom; i++)
    //{
      //  std::cout << "percenent done " << i*1.0/nZ_closeGeom << std::endl;
    typedef std::chrono::high_resolution_clock Time0;
    typedef std::chrono::duration<float> fsec0;
    auto start_clock0 = Time0::now();
       hashGeom geo1(nLines, boundaries.data(), 
                        closeGeomGridr.data(), closeGeomGridy.data(),
                        closeGeomGridz.data(),
                        n_closeGeomElements, closeGeom.data(),
                        nR_closeGeom,nY_closeGeom,nZ_closeGeom);
       thrust::for_each(thrust::device, lines0+world_rank*nHashMeshPointsPerProcess,lines0+world_rank*nHashMeshPointsPerProcess+hashMeshIncrements[world_rank]-1,geo1);
       //for(int i=0;i<nR_closeGeom*nY_closeGeom*nZ_closeGeom;i++)
       //{
       // geo1(i);
       //}
       #if USE_CUDA
         cudaDeviceSynchronize();
       #endif
    //}
#if USE_MPI > 0
    MPI_Barrier(MPI_COMM_WORLD);
    //Collect stuff
   for(int rr=1; rr<world_size;rr++)
{
if(world_rank == rr)
{
    MPI_Send(&closeGeom[world_rank*nHashMeshPointsPerProcess*n_closeGeomElements], hashMeshIncrements[world_rank]*n_closeGeomElements, MPI_INT, 0, 0, MPI_COMM_WORLD);
}
else if(world_rank == 0)
{
    MPI_Recv(&closeGeom[rr*nHashMeshPointsPerProcess*n_closeGeomElements], hashMeshIncrements[rr]*n_closeGeomElements, MPI_INT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(closeGeom.data(), nR_closeGeom*nY_closeGeom*nZ_closeGeom*n_closeGeomElements,MPI_INT,0,MPI_COMM_WORLD);
#endif
    #if USE_CUDA
      cudaDeviceSynchronize();
    #endif
      auto finish_clock0 = Time0::now();
      fsec0 fs0 = finish_clock0 - start_clock0;
      printf("Time taken          is %6.3f (secs) \n", fs0.count());
#if USE_MPI > 0  
      if(world_rank == 0)
      {
#endif
    NcFile ncFile_hash("output/geomHash.nc", NcFile::replace);
    NcDim hashNR = ncFile_hash.addDim("nR",nR_closeGeom);
    #if USE3DTETGEOM > 0
      NcDim hashNY = ncFile_hash.addDim("nY",nY_closeGeom);
    #endif
    NcDim hashNZ = ncFile_hash.addDim("nZ",nZ_closeGeom);
    NcDim hashN = ncFile_hash.addDim("n",n_closeGeomElements);
    vector<NcDim> geomHashDim;
    geomHashDim.push_back(hashNR);
    #if USE3DTETGEOM > 0
      geomHashDim.push_back(hashNY);
    #endif
    geomHashDim.push_back(hashNZ);
    geomHashDim.push_back(hashN);
    NcVar hash_gridR = ncFile_hash.addVar("gridR",ncFloat,hashNR);
    #if USE3DTETGEOM > 0
      NcVar hash_gridY = ncFile_hash.addVar("gridY",ncFloat,hashNY);
    #endif
    NcVar hash_gridZ = ncFile_hash.addVar("gridZ",ncFloat,hashNZ);
    NcVar hash = ncFile_hash.addVar("hash",ncInt,geomHashDim);
    hash_gridR.putVar(&closeGeomGridr[0]);
    #if USE3DTETGEOM > 0
      hash_gridY.putVar(&closeGeomGridy[0]);
    #endif

    hash_gridZ.putVar(&closeGeomGridz[0]);
    hash.putVar(&closeGeom[0]);
    ncFile_hash.close();
#if USE_MPI > 0
      }
#endif
  #elif GEOM_HASH > 1
    getVarFromFile(cfg,input_path+hashFile,geomHashCfg,"gridRString",closeGeomGridr);
    getVarFromFile(cfg,input_path+hashFile,geomHashCfg,"gridZString",closeGeomGridz);
    #if USE3DTETGEOM >0
      getVarFromFile(cfg,input_path+hashFile,geomHashCfg,"gridYString",closeGeomGridy);
    #endif
    std::cout << "geom hash numbers " << nR_closeGeom << " " << nY_closeGeom
              << " " << nZ_closeGeom << " " << n_closeGeomElements << std::endl;
    getVarFromFile(cfg,input_path+hashFile,geomHashCfg,"closeGeomString",closeGeom);
    std::cout << "geom hash indices " << closeGeom[0] << " " << closeGeom[1] << " " << closeGeom[2] << std::endl;
  #endif
              
  int nR_closeGeom_sheath = 1;
  int nY_closeGeom_sheath = 1;
  int nZ_closeGeom_sheath = 1;
  int n_closeGeomElements_sheath = 1;
  int nGeomHash_sheath = 1;
  std::string geomHashSheathCfg= "geometry_sheath.";  
  #if GEOM_HASH_SHEATH == 1
    getVariable(cfg,geomHashSheathCfg+"nR_closeGeom",nR_closeGeom_sheath);
    getVariable(cfg,geomHashSheathCfg+"nZ_closeGeom",nZ_closeGeom_sheath);
    getVariable(cfg,geomHashSheathCfg+"n_closeGeomElements",n_closeGeomElements_sheath);
    nGeomHash_sheath = nR_closeGeom_sheath*nZ_closeGeom_sheath*n_closeGeomElements_sheath;
    #if USE3DTETGEOM > 0
      getVariable(cfg,geomHashSheathCfg+"nY_closeGeom",nY_closeGeom_sheath);
      nGeomHash_sheath = nY_closeGeom_sheath*nGeomHash_sheath;
    #endif
  #endif

  #if GEOM_HASH_SHEATH > 1
    std::string hashFile_sheath;
    getVariable(cfg,geomHashSheathCfg+"fileString",hashFile_sheath);
    nR_closeGeom_sheath = getDimFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridNrString");
    nZ_closeGeom_sheath = getDimFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridNzString");
    n_closeGeomElements_sheath = getDimFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"nearestNelementsString");
    nGeomHash_sheath = nR_closeGeom_sheath*nZ_closeGeom_sheath*n_closeGeomElements_sheath;
    #if USE3DTETGEOM > 0
      nY_closeGeom_sheath = getDimFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridNyString");
      nGeomHash_sheath = nY_closeGeom_sheath*nGeomHash_sheath;
    #else
    #endif
  #endif
  sim::Array<float> closeGeomGridr_sheath(nR_closeGeom_sheath), 
                    closeGeomGridy_sheath(nY_closeGeom_sheath),
                    closeGeomGridz_sheath(nZ_closeGeom_sheath);
  sim::Array<int> closeGeom_sheath(nGeomHash_sheath);
  #if GEOM_HASH_SHEATH  ==1
    float hashX0_s,hashX1_s,hashY0_s,hashY1_s,hashZ0_s,hashZ1_s;
    getVariable(cfg,geomHashSheathCfg+"hashX0",hashX0_s);
    getVariable(cfg,geomHashSheathCfg+"hashX1",hashX1_s);
    getVariable(cfg,geomHashSheathCfg+"hashZ0",hashZ0_s);
    getVariable(cfg,geomHashSheathCfg+"hashZ1",hashZ1_s);
    #if USE3DTETGEOM > 0
      getVariable(cfg,geomHashSheathCfg+"hashY0",hashY0_s);
      getVariable(cfg,geomHashSheathCfg+"hashY1",hashY1_s);
    #endif
    
    for(int i=0; i<nR_closeGeom_sheath; i++)
    {  closeGeomGridr_sheath[i] = (hashX1_s - hashX0_s)*i/(nR_closeGeom_sheath - 1)+ hashX0_s;}
    for(int j=0; j<nY_closeGeom_sheath; j++)
    {  closeGeomGridy_sheath[j] = (hashY1_s - hashY0_s)*j/(nY_closeGeom_sheath - 1)+ hashY0_s;}
    for(int k=0; k<nZ_closeGeom_sheath; k++)
    {  closeGeomGridz_sheath[k] = (hashZ1_s - hashZ0_s)*k/(nZ_closeGeom_sheath - 1)+ hashZ0_s;}
      
    thrust::counting_iterator<std::size_t> lines0_s(0);  
    thrust::counting_iterator<std::size_t> lines1_s(nR_closeGeom_sheath*nY_closeGeom_sheath);
    sim::Array<float> minDist1_s(nGeomHash_sheath,1e6);
    int nHashMeshPointsPerProcess_s=ceil(nR_closeGeom_sheath*nY_closeGeom_sheath*nZ_closeGeom_sheath/world_size);
   std::vector<int> hashMeshIncrements_s(world_size);
   for(int j=0;j<world_size-1;j++)
   {
       hashMeshIncrements_s[j] = nHashMeshPointsPerProcess_s;
   }
   hashMeshIncrements_s[world_size-1]=nR_closeGeom_sheath*nY_closeGeom_sheath*nZ_closeGeom_sheath - (world_size-1)*nHashMeshPointsPerProcess_s;
    typedef std::chrono::high_resolution_clock Time0_s;
    typedef std::chrono::duration<float> fsec0_s;
    auto start_clock0_s = Time0_s::now();
       hashGeom_sheath geo_s(nLines, boundaries.data(), 
                        closeGeomGridr_sheath.data(), closeGeomGridy_sheath.data(),
                        closeGeomGridz_sheath.data(),
                        n_closeGeomElements_sheath, closeGeom_sheath.data(),
                        nR_closeGeom_sheath,nY_closeGeom_sheath,nZ_closeGeom_sheath);
       thrust::for_each(thrust::device, lines0_s+world_rank*nHashMeshPointsPerProcess_s,lines0_s+world_rank*nHashMeshPointsPerProcess_s+hashMeshIncrements_s[world_rank]-1,geo_s);
       #if USE_CUDA
         cudaDeviceSynchronize();
       #endif
#if USE_MPI > 0
    MPI_Barrier(MPI_COMM_WORLD);
    //Collect stuff
   for(int rr=1; rr<world_size;rr++)
{
if(world_rank == rr)
{
    MPI_Send(&closeGeom_sheath[world_rank*nHashMeshPointsPerProcess_s*n_closeGeomElements_sheath], hashMeshIncrements_s[world_rank]*n_closeGeomElements_sheath, MPI_INT, 0, 0, MPI_COMM_WORLD);
}
else if(world_rank == 0)
{
    MPI_Recv(&closeGeom_sheath[rr*nHashMeshPointsPerProcess_s*n_closeGeomElements_sheath], hashMeshIncrements_s[rr]*n_closeGeomElements_sheath, MPI_INT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(closeGeom_sheath.data(), nR_closeGeom_sheath*nY_closeGeom_sheath*nZ_closeGeom_sheath*n_closeGeomElements_sheath,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    #if USE_CUDA
      cudaDeviceSynchronize();
    #endif
      auto finish_clock0_s = Time0_s::now();
      fsec0_s fs0_s = finish_clock0_s - start_clock0_s;
      printf("Time taken          is %6.3f (secs) \n", fs0_s.count());
#if USE_MPI > 0  
      if(world_rank == 0)
      {
#endif
    
    NcFile ncFile_hash_sheath("output/geomHash_sheath.nc", NcFile::replace);
    NcDim hashNR_sheath = ncFile_hash_sheath.addDim("nR",nR_closeGeom_sheath);
    NcDim hashNY_sheath = ncFile_hash_sheath.addDim("nY",nY_closeGeom_sheath);
    NcDim hashNZ_sheath = ncFile_hash_sheath.addDim("nZ",nZ_closeGeom_sheath);
    NcDim hashN_sheath = ncFile_hash_sheath.addDim("n",n_closeGeomElements_sheath);
    vector<NcDim> geomHashDim_sheath;
    geomHashDim_sheath.push_back(hashNR_sheath);
    geomHashDim_sheath.push_back(hashNY_sheath);
    geomHashDim_sheath.push_back(hashNZ_sheath);
    geomHashDim_sheath.push_back(hashN_sheath);
    NcVar hash_gridR_sheath = ncFile_hash_sheath.addVar("gridR",ncDouble,hashNR_sheath);
    NcVar hash_gridY_sheath = ncFile_hash_sheath.addVar("gridY",ncDouble,hashNY_sheath);
    NcVar hash_gridZ_sheath = ncFile_hash_sheath.addVar("gridZ",ncDouble,hashNZ_sheath);
    NcVar hash_sheath = ncFile_hash_sheath.addVar("hash",ncDouble,geomHashDim_sheath);
    hash_gridR_sheath.putVar(&closeGeomGridr_sheath[0]);
    hash_gridY_sheath.putVar(&closeGeomGridy_sheath[0]);
    hash_gridZ_sheath.putVar(&closeGeomGridz_sheath[0]);
    hash_sheath.putVar(&closeGeom_sheath[0]);
    ncFile_hash_sheath.close();
#if USE_MPI > 0  
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    #if USE_CUDA
      cudaDeviceSynchronize();
    #endif
  #elif GEOM_HASH_SHEATH > 1
    getVarFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridRString",closeGeomGridr_sheath);
    getVarFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridZString",closeGeomGridz_sheath);
    #if USE3DTETGEOM >0
      getVarFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"gridYString",closeGeomGridy_sheath);
    #endif
      getVarFromFile(cfg,input_path+hashFile_sheath,geomHashSheathCfg,"closeGeomString",closeGeom_sheath);
  #endif

  int nR_Lc = 1;
  int nY_Lc = 1;
  int nZ_Lc = 1;
  int nTracers = 1;
  std::string connLengthCfg= "connectionLength.";  
  std::string lcFile;
  getVariable(cfg,connLengthCfg+"fileString",lcFile);
  #if GENERATE_LC == 1
    getVariable(cfg,connLengthCfg+"nX",nR_Lc);
    getVariable(cfg,connLengthCfg+"nY",nY_Lc);
    getVariable(cfg,connLengthCfg+"nZ",nZ_Lc);
    float r0_Lc, r1_Lc, y0_Lc,y1_Lc,z0_Lc,z1_Lc, dr;
    int nTraceSteps;
    getVariable(cfg,connLengthCfg+"netx0",r0_Lc);
    getVariable(cfg,connLengthCfg+"netx1",r1_Lc);
    getVariable(cfg,connLengthCfg+"nety0",y0_Lc);
    getVariable(cfg,connLengthCfg+"nety1",y1_Lc);
    getVariable(cfg,connLengthCfg+"netz0",z0_Lc);
    getVariable(cfg,connLengthCfg+"netz1",z1_Lc);
    getVariable(cfg,connLengthCfg+"nTraceSteps",nTraceSteps);
    getVariable(cfg,connLengthCfg+"dr",dr);
  #endif 
  #if LC_INTERP > 1
    nR_Lc = getDimFromFile(cfg,input_path+lcFile,connLengthCfg,"gridNrString");
    nY_Lc = getDimFromFile(cfg,input_path+lcFile,connLengthCfg,"gridNyString");
    nZ_Lc = getDimFromFile(cfg,input_path+lcFile,connLengthCfg,"gridNzString");
  #endif
  
  #if USE3DTETGEOM > 0
    nTracers = nR_Lc*nY_Lc*nZ_Lc;
  #else
    nTracers = nR_Lc*nZ_Lc;
  #endif
  
  sim::Array<float> Lc(nTracers),s(nTracers);
  sim::Array<float> gridRLc(nR_Lc),gridYLc(nY_Lc),gridZLc(nZ_Lc);
  sim::Array<int> noIntersectionNodes(nTracers);
  #if GENERATE_LC ==1
  float lcBuffer = 0.75;
   // if( !boost::filesystem::exists( lcFile ) )
   // {
    //   std::cout << "No pre-existing connection length file found" << std::endl;    
      #if USE3DTETGEOM > 0
        float dy_Lc = (y1_Lc-y0_Lc)/(nY_Lc-1);
        for(int j=0;j<nY_Lc; j++)
        {
         gridYLc[j] = y0_Lc + j*dy_Lc;
        }
      #endif
      float dr_Lc = (r1_Lc-r0_Lc)/(nR_Lc-1);
      for(int i=0;i<nR_Lc; i++)
      {
       gridRLc[i] = r0_Lc + i*dr_Lc;
      }

      float dz_Lc = (z1_Lc-z0_Lc)/(nZ_Lc-1);
      for(int j=0;j<nZ_Lc; j++)
      {
       gridZLc[j] = z0_Lc + j*dz_Lc;
      }
      std::cout << "Creating tracer particles" << std::endl;
      thrust::counting_iterator<std::size_t> lcBegin(0);  
      thrust::counting_iterator<std::size_t> lcEnd(nTracers);
      auto forwardTracerParticles = new Particles(nTracers);
      auto backwardTracerParticles = new Particles(nTracers);
      int addIndex = 0;
      std::cout << "Initializing tracer particles" << std::endl;
   
      for(int i=0;i<nR_Lc; i++)
      {
        for(int j=0;j<nY_Lc;j++)
        {
           for (int k=0;k<nZ_Lc;k++)
           {
              #if USE3DTETGEOM > 0
                addIndex = i + j*nR_Lc + k*nR_Lc*nY_Lc;
              #else
                addIndex = i+k*nR_Lc;
              #endif
              forwardTracerParticles->setParticle(addIndex,gridRLc[i], gridYLc[j], gridZLc[k], 
                                                  0.0, 0.0, 0.0, 0, 0.0, 0.0);
              backwardTracerParticles->setParticle(addIndex,gridRLc[i], gridYLc[j], gridZLc[k], 
                                                   0.0, 0.0, 0.0, 0, 0.0, 0.0);
           }
        }
      }
      
      //dummy surfaces for Lcs calculation (geometry_check)
      auto dummy_surfaces = new Surfaces(1,1,1);
      dummy_surfaces->setSurface(1, 1,1,1 ,1 ,1);
                              
      typedef std::chrono::high_resolution_clock Time_trace;
      typedef std::chrono::duration<float> fsec_trace;
      auto start_clock_trace = Time_trace::now();
      std::cout << "Starting trace loop" << std::endl;
      std::cout << "nTraceSteps"<< nTraceSteps << " dr "<< dr  << std::endl;
      for (int ii=0;ii<nTraceSteps; ii++)
      {
        #if USE_CUDA 
          cudaDeviceSynchronize();
        #endif
        thrust::for_each(thrust::device, lcBegin,lcEnd,
                 field_line_trace(1.0,forwardTracerParticles,dr,boundaries.data(), nLines,
                     nR_Lc,nZ_Lc,gridRLc.data(),gridZLc.data(),Lc.data(),
                     nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                     &br.front(),&bz.front(),&by.front()));
        
        thrust::for_each(thrust::device, lcBegin,lcEnd,
                 field_line_trace(-1.0,backwardTracerParticles,dr,boundaries.data(), nLines,
                     nR_Lc,nZ_Lc,gridRLc.data(),gridZLc.data(),Lc.data(),
                     nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                     &br.front(),&bz.front(),&by.front()));
            
        thrust::for_each(thrust::device, lcBegin,lcEnd,
                  geometry_check(forwardTracerParticles,nLines,&boundaries[0], dummy_surfaces,dr,ii,
                      nR_closeGeom,nY_closeGeom,nZ_closeGeom,n_closeGeomElements,
                      &closeGeomGridr.front(),&closeGeomGridy.front(),&closeGeomGridz.front(),
                      &closeGeom.front()) );
        
        thrust::for_each(thrust::device, lcBegin,lcEnd,
                  geometry_check(backwardTracerParticles,nLines,&boundaries[0],dummy_surfaces,dr,ii,
                      nR_closeGeom,nY_closeGeom,nZ_closeGeom,n_closeGeomElements,
                      &closeGeomGridr.front(),&closeGeomGridy.front(),&closeGeomGridz.front(),
                      &closeGeom.front()) );
      }
      auto finish_clock_trace = Time_trace::now();
      fsec_trace fstrace = finish_clock_trace - start_clock_trace;
      printf("Time taken          is %6.3f (secs) \n", fstrace.count());
      printf("Time taken per step is %6.3f (secs) \n", fstrace.count() / (float) nTraceSteps);
      #if USE_CUDA 
         cudaDeviceSynchronize();
      #endif
      addIndex = 0;
      float forwardDist = 0.0;
      float backwardDist = 0.0;
      for(int i=0;i<nR_Lc; i++)
      {
        for(int j=0;j<nY_Lc;j++)
        {
          for(int k=0;k<nZ_Lc;k++)
          {
              //std::cout << "hitwall of tracers " << forwardTracerParticles->hitWall[addIndex] << std::endl; 
             #if USE3DTETGEOM > 0
               addIndex = i + j*nR_Lc + k*nR_Lc*nY_Lc;
             #else
                   addIndex = i+k*nR_Lc;
             #endif
             if(forwardTracerParticles->hitWall[addIndex] > 0)
             {
                forwardDist = forwardTracerParticles->distanceTraveled[addIndex];
             }
             else
             { forwardDist = 0.0;}

             if(backwardTracerParticles->hitWall[addIndex] > 0)
             {
                backwardDist = backwardTracerParticles->distanceTraveled[addIndex];
             }
             else backwardDist = 0.0;
             
             Lc[addIndex] = forwardDist + backwardDist;
             if(Lc[addIndex] > 0.0)
             {
                Lc[addIndex] = Lc[addIndex] + lcBuffer;
             } 
      //       if(forwardTracerParticles->distanceTraveled[addIndex] > 
      //               backwardTracerParticles->distanceTraveled[addIndex])
             if(forwardTracerParticles->distanceTraveled[addIndex] > 
                     0.5*Lc[addIndex])
             {
               s[addIndex] = -(0.5*Lc[addIndex]-backwardTracerParticles->distanceTraveled[addIndex]);
             }
             else
             {
               s[addIndex] = (0.5*Lc[addIndex]-forwardTracerParticles->distanceTraveled[addIndex]);
             }
             if(forwardTracerParticles->hitWall[addIndex] + backwardTracerParticles->hitWall[addIndex]<4.0)
             {
               noIntersectionNodes[addIndex] = 1;    
             }
          }
        }
      }

      NcFile ncFileLC("LcS.nc", NcFile::replace);
      NcDim nc_nTracers = ncFileLC.addDim("nTracers",nTracers);
      NcDim nc_nRLc = ncFileLC.addDim("nR",nR_Lc);
      NcDim nc_nYLc = ncFileLC.addDim("nY",nY_Lc);
      NcDim nc_nZLc = ncFileLC.addDim("nZ",nZ_Lc);
      
      NcVar nc_Lc = ncFileLC.addVar("Lc",ncDouble,nc_nTracers);
      NcVar nc_s = ncFileLC.addVar("s",ncDouble,nc_nTracers);
      NcVar nc_nI = ncFileLC.addVar("noIntersection",ncDouble,nc_nTracers);
      NcVar nc_gridRLc = ncFileLC.addVar("gridR",ncDouble,nc_nRLc);
      NcVar nc_gridYLc = ncFileLC.addVar("gridY",ncDouble,nc_nYLc);
      NcVar nc_gridZLc = ncFileLC.addVar("gridZ",ncDouble,nc_nZLc);
      
      nc_Lc.putVar(&Lc[0]);
      nc_s.putVar(&s[0]);
      nc_nI.putVar(&noIntersectionNodes[0]);
      nc_gridRLc.putVar(&gridRLc[0]);
      nc_gridYLc.putVar(&gridYLc[0]);
      nc_gridZLc.putVar(&gridZLc[0]);
      ncFileLC.close();
  #if USE_CUDA 
             cudaDeviceSynchronize();
      #endif
   //}         
  #endif    

  #if LC_INTERP > 1
    std::cout << "Importing pre-existing connection length file" << std::endl;
    getVariable(cfg,connLengthCfg+"fileString",lcFile);
    getVarFromFile(cfg,input_path+lcFile,connLengthCfg,"gridRString",gridRLc);
    getVarFromFile(cfg,input_path+lcFile,connLengthCfg,"gridYString",gridYLc);
    getVarFromFile(cfg,input_path+lcFile,connLengthCfg,"gridZString",gridZLc);
    getVarFromFile(cfg,input_path+lcFile,connLengthCfg,"LcString",Lc);
    getVarFromFile(cfg,input_path+lcFile,connLengthCfg,"SString",s);
  #endif
  
  //Background Plasma Temperature Initialization    
  int nR_Temp = 1;
  int nY_Temp = 1;
  int nZ_Temp = 1;
  int n_Temp = 1;
  std::string tempCfg = "backgroundPlasmaProfiles.Temperature.";
  #if TEMP_INTERP > 0
    std::string tempFile;
    getVariable(cfg,tempCfg+"fileString",tempFile);
    nR_Temp = getDimFromFile(cfg,input_path+tempFile,tempCfg,"gridNrString");
  #endif
  #if TEMP_INTERP > 1
    nZ_Temp = getDimFromFile(cfg,input_path+tempFile,tempCfg,"gridNzString");
  #endif
  #if TEMP_INTERP > 2
    nY_Temp = getDimFromFile(cfg,input_path+tempFile,tempCfg,"gridNyString");
  #endif
  
  sim::Array<float> TempGridr(nR_Temp), TempGridz(nZ_Temp), TempGridy(nY_Temp);
  n_Temp = nR_Temp*nY_Temp*nZ_Temp;
  sim::Array<float> ti(n_Temp), te(n_Temp);

  #if TEMP_INTERP == 0
    getVariable(cfg,tempCfg+"ti",ti[0]);
    getVariable(cfg,tempCfg+"te",te[0]);
  #else
    getVarFromFile(cfg,input_path+tempFile,tempCfg,"gridRString",TempGridr);
    #if TEMP_INTERP > 1
      getVarFromFile(cfg,input_path+tempFile,tempCfg,"gridZString",TempGridz);
    #endif
    #if TEMP_INTERP > 2
      getVarFromFile(cfg,input_path+tempFile,tempCfg,"gridYString",TempGridy);
    #endif
    getVarFromFile(cfg,input_path+tempFile,tempCfg,"IonTempString",ti);
    getVarFromFile(cfg,input_path+tempFile,tempCfg,"ElectronTempString",te);
  #endif

  float testVec = 0.0;
  testVec = interp2dCombined(0.01,0.0,0.0,nR_Temp,
                    nZ_Temp,TempGridr.data(),TempGridz.data(),ti.data());
  std::cout << "Finished Temperature import "<< testVec << std::endl; 
  
  //Background Plasma Density Initialization    
  int nR_Dens = 1;
  int nY_Dens = 1;
  int nZ_Dens = 1;
  int n_Dens = 1;
  std::string densCfg = "backgroundPlasmaProfiles.Density.";
  #if DENSITY_INTERP > 0
    std::string densFile;
    getVariable(cfg,densCfg+"fileString",densFile);
    nR_Dens = getDimFromFile(cfg,input_path+densFile,densCfg,"gridNrString");
  #endif
  #if DENSITY_INTERP > 1
    nZ_Dens = getDimFromFile(cfg,input_path+densFile,densCfg,"gridNzString");
  #endif
  #if DENSITY_INTERP > 2
    nY_Dens = getDimFromFile(cfg,input_path+densFile,densCfg,"gridNyString");
  #endif
  
  sim::Array<float> DensGridr(nR_Dens), DensGridz(nZ_Dens), DensGridy(nY_Dens);
  n_Dens = nR_Dens*nY_Dens*nZ_Dens;
  sim::Array<float> ni(n_Dens), ne(n_Dens);

  #if DENSITY_INTERP == 0
    getVariable(cfg,densCfg+"ni",ni[0]);
    getVariable(cfg,densCfg+"ne",ne[0]);
  #else
    getVarFromFile(cfg,input_path+densFile,densCfg,"gridRString",DensGridr);
    #if DENSITY_INTERP > 1
      getVarFromFile(cfg,input_path+densFile,densCfg,"gridZString",DensGridz);
    #endif
    #if DENSITY_INTERP > 2
      getVarFromFile(cfg,input_path+densFile,densCfg,"gridYString",DensGridy);
    #endif
    getVarFromFile(cfg,input_path+densFile,densCfg,"IonDensityString",ni);
    getVarFromFile(cfg,input_path+densFile,densCfg,"ElectronDensityString",ne);
  #endif
  std::cout << "Finished density import "<< interp2dCombined(0.001,0.0,0.1,nR_Dens,nZ_Dens,
                         &DensGridr.front(),&DensGridz.front(),&ne.front())
 <<" " << interp2dCombined(0.02,0.0,0.1,nR_Dens,nZ_Dens,
                                 &DensGridr.front(),&DensGridz.front(),&ne.front()) << std::endl;
  for(int i=0;i<100;i++)
  {
      std::cout << i*0.001 << " " << interp2dCombined(0.001*i,0.0,0.0,nR_Dens,nZ_Dens,
                                       &DensGridr.front(),&DensGridz.front(),&ne.front()) << std::endl;
  }
  std::cout << " z=0.1" << std::endl;
  for(int i=0; i<100;i++)
  {
      std::cout << i*0.001 << " " << interp2dCombined(0.001*i,0.0,0.1,nR_Dens,nZ_Dens,
                                       &DensGridr.front(),&DensGridz.front(),&ne.front()) << std::endl;
  }
  //Background Plasma flow velocity initialization    
  int nR_flowV = 1;
  int nY_flowV = 1;
  int nZ_flowV = 1;
  int n_flowV = 1;
  std::string flowVCfg = "backgroundPlasmaProfiles.FlowVelocity.";
  #if FLOWV_INTERP == 1
    nR_flowV=nR_Lc;
    nY_flowV=nY_Lc;
    nZ_flowV=nZ_Lc;
  #endif
  #if FLOWV_INTERP > 1
    std::string flowVFile;
    getVariable(cfg,flowVCfg+"fileString",flowVFile);
    nR_flowV = getDimFromFile(cfg,input_path+flowVFile,flowVCfg,"gridNrString");
    nZ_flowV = getDimFromFile(cfg,input_path+flowVFile,flowVCfg,"gridNzString");
  #endif
  #if FLOWV_INTERP > 2
    nY_flowV = getDimFromFile(cfg,input_path+flowVFile,flowVCfg,"gridNyString");
  #endif

  sim::Array<float> flowVGridr(nR_flowV),flowVGridy(nY_flowV), flowVGridz(nZ_flowV);
  n_flowV = nR_flowV*nY_flowV*nZ_flowV;
  sim::Array<float> flowVr(n_flowV), flowVz(n_flowV),flowVt(n_flowV);

  #if FLOWV_INTERP == 0
    getVariable(cfg,flowVCfg+"flowVr",flowVr[0]);
    getVariable(cfg,flowVCfg+"flowVy",flowVt[0]);
    getVariable(cfg,flowVCfg+"flowVz",flowVz[0]);
  #else
    #if FLOWV_INTERP > 1
      getVarFromFile(cfg,input_path+flowVFile,flowVCfg,"flowVrString",flowVr);
      getVarFromFile(cfg,input_path+flowVFile,flowVCfg,"flowVtString",flowVt);
    #endif
    #if FLOWV_INTERP > 2  
      getVarFromFile(cfg,input_path+flowVFile,flowVCfg,"flowVzString",flowVz);
    #endif
  #endif
  #if FLOWV_INTERP == 1
    for(int i=0;i<nR_flowV;i++)
    {
      std::cout << " !!! gridRLc " << gridRLc[i] << std::endl;
    }
    std::cout << " !!! gridZLc " << gridZLc[0] << std::endl;
    for(int i=0;i<nR_flowV;i++)
    {
        flowVGridr[i] = gridRLc[i];
    }
    for(int i=0;i<nZ_flowV;i++)
    {
        flowVGridz[i] = gridZLc[i];
    }
    std::cout << " !!! flowvgridr0 " << flowVGridr[0] << std::endl;
    int nFlowVs = nR_Lc*nZ_Lc;
    #if LC_INTERP == 3
      for(int i=0;i<nY_flowV;i++)
        flowVGridy[i] = gridYLc[i];
      nFlowVs = nR_Lc*nY_Lc*nZ_Lc;
    #endif
    float thisY = 0.0;
    float cs0=0.0;
    float teLocal = 0.0;
    float tiLocal = 0.0;
    float BLocal[3] = {0.0,0.0,0.0};
    float Bnorm[3] = {0.0,0.0,0.0};
    float Bmag = 0.0;
    int index = 0;
    float cs = 0.0;
    float absS = 0.0;
    std::cout << "Beginning analytic flowV calculation "<< std::endl; 
    for(int i=0;i<nR_Lc; i++)
    {
     #if LC_INTERP == 3
     for(int k=0;k < nY_Lc; k++)
     { thisY = flowVGridy[k];
     #endif

         for(int j=0;j<nZ_Lc;j++)
      { 
        //std::cout << "debug here 1 " << i << " " << j << std::endl;
        //std::cout << "debug here 2 " << flowVGridr[i] << " " << thisY << " "  
        //    << flowVGridz[j] << " " << nR_Temp << " "<<nZ_Temp << std::endl;
        teLocal = interp2dCombined(flowVGridr[i],thisY,flowVGridz[j],nR_Temp,nZ_Temp, 
                &TempGridr.front(),&TempGridz.front(),&te.front());
        tiLocal = interp2dCombined(flowVGridr[i],thisY,flowVGridz[j],nR_Temp,nZ_Temp, 
                &TempGridr.front(),&TempGridz.front(),&ti.front());
        cs0 = sqrt((teLocal+tiLocal)*1.602e-19/(background_amu*1.66e-27));
        interp2dVector(&BLocal[0],flowVGridr[i],thisY,flowVGridz[j],nR_Bfield,
                    nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
        Bmag = sqrt(BLocal[0]*BLocal[0] + BLocal[1]*BLocal[1] + BLocal[2]*BLocal[2]);
        Bnorm[0] = BLocal[0]/Bmag;
        Bnorm[1] = BLocal[1]/Bmag;
        Bnorm[2] = BLocal[2]/Bmag;

     #if LC_INTERP == 3
        index = i+k*nR_Lc + j*nR_Lc*nY_Lc;
        //std::cout << "flowv calc index " << index << std::endl;
     #else
        index = i+j*nR_Lc;
     #endif
        absS = abs(s[index]);
        cs = cs0*(0.5*Lc[index]/absS - sqrt(0.25*Lc[index]*Lc[index]/absS/absS - 1.0));
        if(std::isnan(cs)) cs = 0.0;
        flowVr[index] = sgn(s[index])*Bnorm[0]*cs;
        flowVt[index] = sgn(s[index])*Bnorm[1]*cs;
        flowVz[index] = sgn(s[index])*Bnorm[2]*cs;
     #if LC_INTERP == 3
      }
     #endif
      }
    }
    std::cout << "Done with initial calculation, beginning sorting" << std::endl;
    sim::Array<float> flowVrSub(nFlowVs), flowVzSub(nFlowVs),
                        flowVySub(nFlowVs);
    sim::Array<int> noIntersectionNearestMax(nFlowVs);
    float surroundingMinimumR = 0.0;
    float surroundingMinimumY = 0.0;
    float surroundingMinimumZ = 0.0;
    int iterIndex = 0;
    for(int i=0; i<nR_Lc;i++)
    {
        std::cout << "i of " << i << " " << nR_Lc << std::endl;
        for(int j=0;j<nY_Lc;j++)
        {
            for(int k=0;k<nZ_Lc;k++)
            {
               index = i+j*nR_Lc + k*nR_Lc*nY_Lc;
               if(noIntersectionNodes[index] ==1)
               {
                   surroundingMinimumR = 0.0;
                   surroundingMinimumY = 0.0;
                   surroundingMinimumZ = 0.0;
                       for(int ii=i-1; ii<i+2;ii++)
                       {
                         for(int jj=j-1;jj<j+2;jj++)
                         {
                           for(int kk=k-1;kk<k+2;kk++)
                           {
                               iterIndex = ii+jj*nR_Lc + kk*nR_Lc*nY_Lc;
                               if(iterIndex > 0 && iterIndex < nFlowVs)
                               {
                               if(noIntersectionNodes[iterIndex] ==0)
                               {
                                 if(abs(flowVr[iterIndex])>abs(surroundingMinimumR))
                                 {
                                   surroundingMinimumR = flowVr[iterIndex];
                                 }
                                 if(abs(flowVt[iterIndex])>abs(surroundingMinimumY))
                                 {
                                   surroundingMinimumY = flowVt[iterIndex];
                                 }
                                 if(abs(flowVz[iterIndex])>abs(surroundingMinimumZ))
                                 {
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
    for(int i=0;i<nFlowVs;i++)
    {
            if(i ==282839)
            {
                std::cout<< " noIntersectionNodes " << noIntersectionNodes[i] << std::endl;
            }
               if(noIntersectionNodes[i] ==1)
               {
                  flowVr[i] =  flowVrSub[i];
                  flowVt[i] =  flowVySub[i];
                  flowVz[i] =  flowVzSub[i];
               }

    }
    NcFile ncFileFlow("flowV.nc", NcFile::replace);
    NcDim nFlowV = ncFileFlow.addDim("n_flowV",n_flowV);
    NcDim nc_nRflow = ncFileFlow.addDim("nR",nR_flowV);
    NcDim nc_nYflow = ncFileFlow.addDim("nY",nY_flowV);
    NcDim nc_nZflow = ncFileFlow.addDim("nZ",nZ_flowV);
    vector<NcDim> dimsFlowV;
    dimsFlowV.push_back(nc_nZflow);
    dimsFlowV.push_back(nc_nYflow);
    dimsFlowV.push_back(nc_nRflow);
    NcVar nc_flowVr = ncFileFlow.addVar("flowVr",ncDouble,dimsFlowV);
    NcVar nc_flowVt = ncFileFlow.addVar("flowVt",ncDouble,dimsFlowV);
    NcVar nc_flowVz = ncFileFlow.addVar("flowVz",ncDouble,dimsFlowV);
    nc_flowVr.putVar(&flowVr[0]);
    nc_flowVt.putVar(&flowVt[0]);
    nc_flowVz.putVar(&flowVz[0]);
    ncFileFlow.close();
    std::string outnameFlowVr = "flowVr.m";
    std::string outnameFlowVz = "flowVz.m";
    std::string outnameFlowVt = "flowVt.m";
    #if LC_INTERP == 3
      OUTPUT3d(profiles_folder,outnameFlowVr, nR_flowV,nY_flowV, nZ_flowV, &flowVr.front());
      OUTPUT3d(profiles_folder,outnameFlowVz, nR_flowV,nY_flowV, nZ_flowV, &flowVz.front());
      OUTPUT3d(profiles_folder,outnameFlowVt, nR_flowV,nY_flowV, nZ_flowV, &flowVt.front());
    #else
      OUTPUT2d(profiles_folder,outnameFlowVr, nR_flowV, nZ_flowV, &flowVr.front());
      OUTPUT2d(profiles_folder,outnameFlowVz, nR_flowV, nZ_flowV, &flowVz.front());
      OUTPUT2d(profiles_folder,outnameFlowVt, nR_flowV, nZ_flowV, &flowVt.front());
    #endif
  #endif

  //Background plasma temperature gradient field intitialization    
  int nR_gradT = 1;
  int nY_gradT = 1;
  int nZ_gradT = 1;
  int n_gradT = 1;
  std::string gradTCfg = "backgroundPlasmaProfiles.gradT.";
  #if GRADT_INTERP > 0
    std::string gradTFile;
    getVariable(cfg,gradTCfg+"fileString",gradTFile);
    nR_gradT = getDimFromFile(cfg,input_path+gradTFile,gradTCfg,"gridNrString");
  #endif
  #if GRADT_INTERP > 1
    nZ_gradT = getDimFromFile(cfg,input_path+gradTFile,gradTCfg,"gridNzString");
  #endif
  #if GRADT_INTERP > 2
    nY_gradT = getDimFromFile(cfg,input_path+gradTFile,gradTCfg,"gridNyString");
  #endif
  sim::Array<float> gradTGridr(nR_gradT),gradTGridy(nY_gradT), gradTGridz(nZ_gradT);
  n_gradT = nR_gradT*nY_gradT*nZ_gradT;
  sim::Array<float> gradTeR(n_gradT), gradTeZ(n_gradT),
                    gradTeY(n_gradT), gradTiR(n_gradT), 
                    gradTiZ(n_gradT), gradTiY(n_gradT);    

  #if GRADT_INTERP == 0
    getVariable(cfg,gradTCfg+"gradTeR",gradTeR[0]);  
    getVariable(cfg,gradTCfg+"gradTeY",gradTeY[0]);  
    getVariable(cfg,gradTCfg+"gradTeZ",gradTeZ[0]);  
    getVariable(cfg,gradTCfg+"gradTiR",gradTiR[0]);  
    getVariable(cfg,gradTCfg+"gradTiY",gradTiY[0]);  
    getVariable(cfg,gradTCfg+"gradTiZ",gradTiZ[0]);  
  #else
    getVarFromFile(cfg,input_path+gradTFile,gradTCfg,"gridRString",gradTGridr);
    #if GRADT_INTERP > 1
      getVarFromFile(cfg,input_path+gradTFile,gradTCfg,"gridZString",gradTGridz);
    #endif
    #if GRADT_INTERP > 2
      getVarFromFile(cfg,input_path+gradTFile,gradTCfg,"gridYString",gradTGridy);
      getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTeYString",gradTeY); 
      getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTiYString",gradTiY); 
    #endif

    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTiRString",gradTiR); 
    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTiZString",gradTiZ); 
    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTeRString",gradTeR); 
    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTeZString",gradTeZ);
    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTeYString",gradTeY); 
    getVarFromFile(cfg,input_path+gradTFile, gradTCfg,"gradTiYString",gradTiY); 
  #endif 
  float gradTi[3] = {0.0};
  interp2dVector(&gradTi[0],1.45,0.0,-1.2,nR_gradT,nZ_gradT,
                            gradTGridr.data() ,gradTGridz.data() ,gradTiR.data(),
                            gradTiZ.data(), gradTiY.data() ); 
  std::cout << "thermal gradient interpolation gradTi " << gradTi[0] << " " << 
     gradTi[1] << " " << gradTi[2] << " " << std::endl; 

  auto particleExample = new Particles(1);
  particleExample->setParticle(0,1.5,0.0, 0.0, 4.4, 0.0,0.0, 74.0, 184.0, 1.0);   
  thermalForce tf(particleExample,1.0e-7,background_amu,
                    nR_gradT,nZ_gradT,gradTGridr.data(),gradTGridz.data(),
                    gradTiR.data(),gradTiZ.data(), gradTiY.data(), 
                    gradTeR.data(), gradTeZ.data(), gradTeY.data(), 
                    nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&by.front());
  tf(0);
  std::cout << "thermal force values " <<tf.nR_gradT <<" " <<  tf.dv_ITG[0] << std::endl;  
  //Initialization of ionization and recombination coefficients    
  int nCS_Ionize, nCS_Recombine;
  const char *ionizeNcs,*ionizeNDens,*ionizeNTemp,
             *ionizeDensGrid,*ionizeTempGrid,*ionizeRCvarChar,
             *recombNcs,*recombNDens,*recombNTemp,
             *recombDensGrid,*recombTempGrid,*recombRCvarChar;
  std::string ionizeFile,recombFile;
  if(cfg.lookupValue("impurityParticleSource.ionization.fileString", ionizeFile) &&
     cfg.lookupValue("impurityParticleSource.ionization.nChargeStateString",ionizeNcs) &&
     cfg.lookupValue("impurityParticleSource.ionization.DensGridString",ionizeNDens) &&
     cfg.lookupValue("impurityParticleSource.ionization.TempGridString",ionizeNTemp) &&
     cfg.lookupValue("impurityParticleSource.ionization.TempGridVarName",ionizeTempGrid) &&
     cfg.lookupValue("impurityParticleSource.ionization.DensGridVarName",ionizeDensGrid) &&
     cfg.lookupValue("impurityParticleSource.ionization.CoeffVarName",ionizeRCvarChar))
  { std::cout << "Ionization rate coefficient file: " << ionizeFile << std::endl;}
  else
  { std::cout << "ERROR: Could not get ionization string info from input file " << std::endl;}
  if(cfg.lookupValue("impurityParticleSource.recombination.fileString", recombFile) &&
     cfg.lookupValue("impurityParticleSource.recombination.nChargeStateString",recombNcs) &&
     cfg.lookupValue("impurityParticleSource.recombination.DensGridString",recombNDens) &&
     cfg.lookupValue("impurityParticleSource.recombination.TempGridString",recombNTemp) &&
     cfg.lookupValue("impurityParticleSource.recombination.TempGridVarName",recombTempGrid) &&
     cfg.lookupValue("impurityParticleSource.recombination.DensGridVarName",recombDensGrid) &&
     cfg.lookupValue("impurityParticleSource.recombination.CoeffVarName",recombRCvarChar))
  { std::cout << "Recombination rate coefficient file: " << recombFile << std::endl;}
  else
  { std::cout << "ERROR: Could not get ionization string info from input file " << std::endl;}
  int i0 = read_profileNs(input_path+ionizeFile,ionizeNcs,recombNcs,nCS_Ionize, nCS_Recombine);

  int nTemperaturesIonize, nDensitiesIonize;
  int i1 = read_profileNs(input_path+ionizeFile,ionizeNDens,ionizeNTemp,nDensitiesIonize,nTemperaturesIonize);

  sim::Array<float> rateCoeff_Ionization(nCS_Ionize*nTemperaturesIonize*nDensitiesIonize);
  sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
                        gridDensity_Ionization(nDensitiesIonize);

  int i2 = read_profiles(input_path+ionizeFile,nTemperaturesIonize,nDensitiesIonize,ionizeTempGrid, 
                         gridTemperature_Ionization,ionizeDensGrid,gridDensity_Ionization,
                         ionizeRCvarChar,rateCoeff_Ionization);
   
  int nTemperaturesRecombine, nDensitiesRecombine;
  int i3 = read_profileNs(input_path+recombFile,recombNDens,recombNTemp,
                          nDensitiesRecombine,nTemperaturesRecombine);

  sim::Array<float> rateCoeff_Recombination(nCS_Recombine*nTemperaturesRecombine*nDensitiesRecombine);
  sim::Array<float> gridTemperature_Recombination(nTemperaturesRecombine),
                    gridDensity_Recombination(nDensitiesRecombine);

  int i4 = read_profiles(input_path+recombFile,nTemperaturesRecombine,nDensitiesRecombine,
             recombTempGrid,gridTemperature_Recombination,recombDensGrid,
             gridDensity_Recombination,
             recombRCvarChar,rateCoeff_Recombination);


  //Applying background values at material boundaries
  std::for_each(boundaries.begin(), boundaries.end()-1,
            boundary_init(background_Z,background_amu,
            nR_Dens,nZ_Dens,DensGridr.data(),DensGridz.data(),ni.data(),
            nR_Bfield,nZ_Bfield,bfieldGridr.data(),
            bfieldGridz.data(),br.data(),bz.data(), by.data(),
            nR_Temp,nZ_Temp,TempGridr.data(),
            TempGridz.data(),ti.data(), biasPotential ));

   std::cout << "Completed Boundary Init " << std::endl;
  
  //Efield
  int nR_PreSheathEfield = 1;
  int nY_PreSheathEfield = 1;
  int nZ_PreSheathEfield = 1;
  int nPSEs = 1;
  std::string PSECfg = "backgroundPlasmaProfiles.Efield.";
  //sim::Array<float> preSheathEGridy(1);
  #if USEPRESHEATHEFIELD > 0    

   std::cout << "Using presheath Efield " << std::endl;
    #if PRESHEATH_INTERP == 1
      nR_PreSheathEfield=nR_Lc;
      nY_PreSheathEfield=nY_Lc;
      nZ_PreSheathEfield=nZ_Lc;
    #endif
    #if PRESHEATH_INTERP > 1
      std::string efieldFile;
      getVariable(cfg,PSECfg+"fileString",efieldFile);
      nR_PreSheathEfield = getDimFromFile(cfg,input_path+efieldFile,PSECfg,"gridNrString");
      nZ_PreSheathEfield = getDimFromFile(cfg,input_path+efieldFile,PSECfg,"gridNzString");
    #endif
    #if PRESHEATH_INTERP > 2
      nY_PreSheathEfield = getDimFromFile(cfg,input_path+efieldFile,PSECfg,"gridNyString");
    #endif
    nPSEs = nR_PreSheathEfield*nY_PreSheathEfield*nZ_PreSheathEfield;
    sim::Array<float> preSheathEGridr(nR_PreSheathEfield),preSheathEGridy(nY_PreSheathEfield),
                      preSheathEGridz(nZ_PreSheathEfield);
    sim::Array<float> PSEr(nPSEs), PSEz(nPSEs),PSEt(nPSEs);
    #if PRESHEATH_INTERP == 0
      getVariable(cfg,PSECfg+"Er",PSEr[0]);
      getVariable(cfg,PSECfg+"Et",PSEt[0]);
      getVariable(cfg,PSECfg+"Ez",PSEz[0]);
    #elif PRESHEATH_INTERP > 1
      getVarFromFile(cfg,input_path+efieldFile,PSECfg,"gridRString",preSheathEGridr);
      getVarFromFile(cfg,input_path+efieldFile,PSECfg,"gridZString",preSheathEGridz);
      #if PRESHEATH_INTERP > 2
        getVarFromFile(cfg,input_path+efieldFile,PSECfg,"gridYString",preSheathEGridy);
      #endif

      getVarFromFile(cfg,input_path+efieldFile,PSECfg,"radialComponentString",PSEr);
      getVarFromFile(cfg,input_path+efieldFile,PSECfg,"toroidalComponentString",PSEt);
      getVarFromFile(cfg,input_path+efieldFile,PSECfg,"axialComponentString",PSEz);
    #endif  

    #if PRESHEATH_INTERP == 1
    
      for(int i=0;i<nR_PreSheathEfield; i++)
      {preSheathEGridr[i]=gridRLc[i];
          std::cout << "gridRLc " << gridRLc[i] << std::endl;}
      for(int i=0;i<nY_PreSheathEfield; i++)
      {preSheathEGridy[i]=gridYLc[i];}
      for(int i=0;i<nZ_PreSheathEfield; i++)
      {preSheathEGridz[i]=gridZLc[i];}
      std::cout << "length of PSE vec " << nPSEs << std::endl;
      float teLocal1 = 0.0;
      float BLocal1[3] = {0.0,0.0,0.0};
      float Bnorm1[3] = {0.0,0.0,0.0};
      float Bmag1 = 0.0;
      int index1 = 0;
      float absS1 = 0.0;
      float Epar = 0.0;
      for(int i=0;i<nR_PreSheathEfield; i++)
      {
       #if LC_INTERP == 3
       for(int k=0;k < nY_PreSheathEfield; k++)
       { thisY = preSheathEGridy[k];
       #endif
        for(int j=0;j<nZ_PreSheathEfield;j++)
        { 
          teLocal1 = interp2dCombined(preSheathEGridr[i],0.0,preSheathEGridz[j],nR_Temp,nZ_Temp, 
                  &TempGridr.front(),&TempGridz.front(),&te.front());
          interp2dVector(&BLocal1[0],gridRLc[i],0.0,gridZLc[j],nR_Bfield,
                      nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
          Bmag1 = sqrt(BLocal1[0]*BLocal1[0] + BLocal1[1]*BLocal1[1] + BLocal1[2]*BLocal1[2]);
          Bnorm1[0] = BLocal1[0]/Bmag1;
          Bnorm1[1] = BLocal1[1]/Bmag1;
          Bnorm1[2] = BLocal1[2]/Bmag1;

       #if LC_INTERP == 3
          index1 = i+k*nR_PreSheathEfield + j*nR_PreSheathEfield*nY_PreSheathEfield;
          //std::cout << "flowv calc index " << index << std::endl;
       #else
          index1 = i+j*nR_PreSheathEfield;
       #endif
          absS1 = abs(s[index1]);
          Epar = teLocal1*(0.5*Lc[index1]/absS1/sqrt(0.25*Lc[index1]*Lc[index1]-absS1*absS1)-1.0/absS1);
          if(std::isnan(Epar)) Epar = 0.0;
          PSEr[index1] = sgn(s[index1])*Bnorm1[0]*Epar;
          PSEt[index1] = sgn(s[index1])*Bnorm1[1]*Epar;
          PSEz[index1] = sgn(s[index1])*Bnorm1[2]*Epar;
        }
       #if LC_INTERP == 3
       }     
       #endif
      }
      sim::Array<float> PSErSub(nPSEs), PSEzSub(nPSEs),
                          PSEySub(nPSEs);

      for(int i=0; i<nR_Lc;i++)
      {
          for(int j=0;j<nY_Lc;j++)
          {
              for(int k=0;k<nZ_Lc;k++)
              {
                 index = i+j*nR_Lc + k*nR_Lc*nY_Lc;
                 if(noIntersectionNodes[index] ==1)
                 {
                     surroundingMinimumR = 0.0;
                     surroundingMinimumY = 0.0;
                     surroundingMinimumZ = 0.0;
                         for(int ii=i-1; ii<i+2;ii++)
                         {
                           for(int jj=j-1;jj<j+2;jj++)
                           {
                             for(int kk=k-1;kk<k+2;kk++)
                             {
                                 iterIndex = ii+jj*nR_Lc + kk*nR_Lc*nY_Lc;
                                 if(iterIndex > 0 && iterIndex < nFlowVs)
                                 {
                                 if(noIntersectionNodes[iterIndex] ==0)
                                 {
                                   if(abs(PSEr[iterIndex])>abs(surroundingMinimumR))
                                   {
                                     surroundingMinimumR = PSEr[iterIndex];
                                   }
                                   if(abs(PSEt[iterIndex])>abs(surroundingMinimumY))
                                   {
                                     surroundingMinimumY = PSEt[iterIndex];
                                   }
                                   if(abs(PSEz[iterIndex])>abs(surroundingMinimumZ))
                                   {
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
      for(int i=0;i<nPSEs;i++)
      {
              if(i ==282839)
              {
                  std::cout<< " noIntersectionNodes " << noIntersectionNodes[i] << std::endl;
              }
                 if(noIntersectionNodes[i] ==1)
                 {
                    PSEr[i] =  PSErSub[i];
                    PSEt[i] =  PSEySub[i];
                    PSEz[i] =  PSEzSub[i];
                 }

      }
      NcVar nc_PSEr = ncFileLC.addVar("PSEr",ncDouble,nc_nTracers);
      NcVar nc_PSEt = ncFileLC.addVar("PSEt",ncDouble,nc_nTracers);
      NcVar nc_PSEz = ncFileLC.addVar("PSEz",ncDouble,nc_nTracers);
      nc_PSEr.putVar(&PSEr[0]);
      nc_PSEt.putVar(&PSEt[0]);
      nc_PSEz.putVar(&PSEz[0]);
    #endif
  #else
    nPSEs = nR_PreSheathEfield*nY_PreSheathEfield*nZ_PreSheathEfield;
    sim::Array<float> preSheathEGridr(nR_PreSheathEfield),preSheathEGridy(nY_PreSheathEfield),
                      preSheathEGridz(nZ_PreSheathEfield);
    sim::Array<float> PSEr(nPSEs), PSEz(nPSEs),PSEt(nPSEs);

  #endif
  std::string outnamePSEfieldR = "PSEfieldR.m";
  std::string outnamePSEfieldZ = "PSEfieldZ.m";
  std::string outnamePSEGridR = "PSEgridR.m";
  std::string outnamePSEGridZ = "PSEgridZ.m";
  OUTPUT1d(profiles_folder,outnamePSEGridR, nR_PreSheathEfield, &preSheathEGridr.front());
  OUTPUT1d(profiles_folder,outnamePSEGridZ, nZ_PreSheathEfield, &preSheathEGridz.front());
  
  OUTPUT3d(profiles_folder,outnamePSEfieldR, nR_PreSheathEfield,nY_PreSheathEfield, nZ_PreSheathEfield, &PSEr.front());
  OUTPUT3d(profiles_folder,outnamePSEfieldZ, nR_PreSheathEfield,nY_PreSheathEfield, nZ_PreSheathEfield, &PSEz.front()); 
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
                  thisE, boundaries.data(),nLines,closestBoundaryIndex );
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

#if SPECTROSCOPY > 0
    float netX0=0.0,netX1=0.0,netY0=0.0,netY1=0.0,netZ0=0.0,netZ1=0.0;
    int net_nX=0,net_nY=0,net_nZ=0;
    int nBins=0;

    if(cfg.lookupValue("diagnostics.netx0", netX0) && 
       cfg.lookupValue("diagnostics.netx1", netX1) && 
       cfg.lookupValue("diagnostics.nety0", netY0) && 
       cfg.lookupValue("diagnostics.nety1", netY1) && 
       cfg.lookupValue("diagnostics.netz0", netZ0) && 
       cfg.lookupValue("diagnostics.netz1", netZ1) && 
       cfg.lookupValue("diagnostics.nX", net_nX) && 
       cfg.lookupValue("diagnostics.nY", net_nY) && 
       cfg.lookupValue("diagnostics.nZ", net_nZ) && 
       cfg.lookupValue("diagnostics.densityChargeBins", nBins))
       {std::cout << "Spectroscopy net imported" << std::endl;}
    else
    { std::cout << "ERROR: Could not get spectroscopy net string info from input file " << std::endl;}
    std::cout << "spec bin Ns " << net_nX << " " << net_nY << " " << net_nZ << std::endl; 
    #if SPECTROSCOPY < 3

      sim::Array<float> net_Bins((nBins+1)*net_nX*net_nZ);
      sim::Array<float> net_BinsTotal((nBins+1)*net_nX*net_nZ);
    #else
      sim::Array<float> net_Bins((nBins+1)*net_nX*net_nY*net_nZ);
      sim::Array<float> net_BinsTotal((nBins+1)*net_nX*net_nY*net_nZ);
    #endif

      /*
      for (int i=0; i<nBins*net_nX*net_nZ; i++)
          {
              std::cout << "i " << i << std::endl;
            net_Bins[i] = 0;
              std::cout << "net bins " << net_Bins[i] << std::endl;
            
          }
      */
      sim::Array<float> gridX_bins(net_nX),gridY_bins(net_nY),gridZ_bins(net_nZ);

      for (int i=0; i< net_nX ; i++)
      {
         gridX_bins[i] = netX0 + 1.0/(net_nX-1)*i*(netX1-netX0);
      }
      for (int i=0; i< net_nY ; i++)
      {
         gridY_bins[i] = netY0 + 1.0/(net_nY-1)*i*(netY1-netY0);
      }

      for (int i=0; i< net_nZ ; i++)
      {
         gridZ_bins[i] = netZ0 + i*1.0/(net_nZ-1)*(netZ1-netZ0);
      }
  #endif    

  // Perp DiffusionCoeff initialization - only used when Diffusion interpolator is = 0
  float perpDiffusionCoeff;
  if (cfg.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",perpDiffusionCoeff))
  {}
  else
  {std::cout << "ERROR: could not get perpendicular diffusion coefficient from input file" << std::endl;}
  //Surface model import
  int nE_sputtRefCoeff = 1, nA_sputtRefCoeff=1;
  int nE_sputtRefDistIn = 1, nA_sputtRefDistIn=1;
  int nE_sputtRefDistOut = 1, nA_sputtRefDistOut = 1;
  int nDistE_surfaceModel = 1,nDistA_surfaceModel = 1;
#if USESURFACEMODEL > 0
  std::string surfaceModelCfg = "surfaceModel.";
  std::string surfaceModelFile;
  getVariable(cfg,surfaceModelCfg+"fileString",surfaceModelFile);
  nE_sputtRefCoeff = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefCoeffString");
  nA_sputtRefCoeff = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefCoeffString");
  nE_sputtRefDistIn = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefDistInString");
  nA_sputtRefDistIn = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefDistInString");
  nE_sputtRefDistOut = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nEsputtRefDistOutString");
  nA_sputtRefDistOut = getDimFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"nAsputtRefDistOutString");
  nDistE_surfaceModel = nE_sputtRefDistIn*nA_sputtRefDistIn*nE_sputtRefDistOut;
  nDistA_surfaceModel = nE_sputtRefDistIn*nA_sputtRefDistIn*nA_sputtRefDistOut;
  std::cout <<  " got dimensions of surface model " << std::endl;
#endif
  sim::Array<float> E_sputtRefCoeff(nE_sputtRefCoeff), A_sputtRefCoeff(nA_sputtRefCoeff),
                    Elog_sputtRefCoeff(nE_sputtRefCoeff),
                    energyDistGrid01(nE_sputtRefDistOut),
                    angleDistGrid01(nA_sputtRefDistOut),
                    spyl_surfaceModel(nE_sputtRefCoeff*nA_sputtRefCoeff),
                    rfyl_surfaceModel(nE_sputtRefCoeff*nA_sputtRefCoeff),
                    E_sputtRefDistIn(nE_sputtRefDistIn), A_sputtRefDistIn(nA_sputtRefDistIn),
                    Elog_sputtRefDistIn(nE_sputtRefDistIn),
                    E_sputtRefDistOut(nE_sputtRefDistOut), A_sputtRefDistOut(nA_sputtRefDistOut),
                    ADist_Y(nDistA_surfaceModel),EDist_Y(nDistE_surfaceModel),
                    ADist_R(nDistA_surfaceModel),EDist_R(nDistE_surfaceModel),
                    ADist_CDF_Y(nDistA_surfaceModel),EDist_CDF_Y(nDistE_surfaceModel),
                    ADist_CDF_R(nDistA_surfaceModel),EDist_CDF_R(nDistE_surfaceModel),
                    ADist_CDF_Y_regrid(nDistA_surfaceModel),EDist_CDF_Y_regrid(nDistE_surfaceModel),
                    ADist_CDF_R_regrid(nDistA_surfaceModel),EDist_CDF_R_regrid(nDistE_surfaceModel);
#if USESURFACEMODEL > 0
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefCoeff",E_sputtRefCoeff);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"A_sputtRefCoeff",A_sputtRefCoeff);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefDistIn",E_sputtRefDistIn);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"A_sputtRefDistIn",A_sputtRefDistIn);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"E_sputtRefDistOut",E_sputtRefDistOut);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"A_sputtRefDistOut",A_sputtRefDistOut);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"sputtYldString",spyl_surfaceModel);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"reflYldString",rfyl_surfaceModel);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"EDist_Y",EDist_Y);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"ADist_Y",ADist_Y);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"EDist_R",EDist_R);
  getVarFromFile(cfg,input_path+surfaceModelFile,surfaceModelCfg,"ADist_R",ADist_R);
  //for(int i=0;i<nDistE_surfaceModel;i++)
  //{
  //    std::cout << " Edist diff Y " << EDist_Y[i] << " " << EDist_R[i] << std::endl;
  //}
  for(int i=0;i<nE_sputtRefCoeff;i++)
  {
      Elog_sputtRefCoeff[i] = log10(E_sputtRefCoeff[i]);
  }
  for(int i=0;i<nE_sputtRefDistIn;i++)
  {
      Elog_sputtRefDistIn[i] = log10(E_sputtRefDistIn[i]);
  }
  for(int i=0;i<nE_sputtRefDistOut;i++)
  {
      energyDistGrid01[i] = i*1.0/nE_sputtRefDistOut;
  }
  for(int i=0;i<nA_sputtRefDistOut;i++)
  {
     angleDistGrid01[i] = i*1.0/nA_sputtRefDistOut;
     //std::cout << " angleDistGrid01[i] " << angleDistGrid01[i] << std::endl;
  }
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,EDist_Y.data(),EDist_CDF_Y.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,ADist_Y.data(),ADist_CDF_Y.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,EDist_R.data(),EDist_CDF_R.data());
  make2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,ADist_R.data(),ADist_CDF_R.data());
  for(int k=0;k<nE_sputtRefDistOut;k++)
  {
        std::cout << "Edist_CDF_Y " << EDist_CDF_Y[44*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
  //      std::cout << "cosDist_CDFR " << EDist_CDF_R[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
  }
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,A_sputtRefDistOut[nA_sputtRefDistOut-1],ADist_CDF_Y.data(),ADist_CDF_Y_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,energyDistGrid01.data(),nE_sputtRefDistOut,E_sputtRefDistOut[nE_sputtRefDistOut-1],EDist_CDF_Y.data(),EDist_CDF_Y_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nA_sputtRefDistOut,angleDistGrid01.data(),nA_sputtRefDistOut,A_sputtRefDistOut[nA_sputtRefDistOut-1],ADist_CDF_R.data(),ADist_CDF_R_regrid.data());
 regrid2dCDF(nE_sputtRefDistIn,nA_sputtRefDistIn,nE_sputtRefDistOut,energyDistGrid01.data(),nE_sputtRefDistOut,E_sputtRefDistOut[nE_sputtRefDistOut-1],EDist_CDF_R.data(),EDist_CDF_R_regrid.data());
 // regrid2dCDF(nE_surfaceModel,nA_surfaceModel,nEdistBins_surfaceModel,energyDistGrid01.data(),nEdistBins_surfaceModel,100.0,energyDist_CDF.data(),energyDist_CDFregrid.data());
 for(int k=0;k<nE_sputtRefDistOut;k++)
  {
        std::cout << "EDist_CDFregridY " << EDist_CDF_Y_regrid[44*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
 //       std::cout << "cosDist_CDFregridR " << EDist_CDF_R_regrid[0*nA_sputtRefDistIn*nE_sputtRefDistOut + 0*nE_sputtRefDistOut+k] << std::endl;
  }
  float spylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff, nE_sputtRefCoeff,A_sputtRefCoeff.data(),
                              Elog_sputtRefCoeff.data(),spyl_surfaceModel.data());
  float rfylInterpVal = interp2d(5.0,log10(250.0),nA_sputtRefCoeff, nE_sputtRefCoeff,A_sputtRefCoeff.data(),
                              Elog_sputtRefCoeff.data(),rfyl_surfaceModel.data());
  float spylEInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,ADist_CDF_Y_regrid.data() );
 float sputEInterpVal = interp3d ( 0.44,5.0,log10(250.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
              energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,EDist_CDF_Y_regrid.data() );
  float rflAInterpVal = interp3d ( 0.44,5.0,log10(250.0),nA_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
          angleDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,ADist_CDF_R_regrid.data() );
 float rflEInterpVal = interp3d ( 0.44,5.0,log10(250.0),nE_sputtRefDistOut,nA_sputtRefDistIn,nE_sputtRefDistIn,
              energyDistGrid01.data(),A_sputtRefDistIn.data(),Elog_sputtRefDistIn.data() ,EDist_CDF_R_regrid.data() );
  std::cout << "Finished surface model import " <<spylInterpVal << " " <<  spylEInterpVal << " " << sputEInterpVal << " "<< rfylInterpVal<< " " << rflAInterpVal << " " << rflEInterpVal <<  std::endl; 
#endif
  // Particle time stepping control
  int ionization_nDtPerApply  = cfg.lookup("timeStep.ionization_nDtPerApply");
  int collision_nDtPerApply  = cfg.lookup("timeStep.collision_nDtPerApply");

  #ifdef __CUDACC__
    cout<<"Using THRUST"<<endl;
  #else
    cout<<"Not using THRUST"<<endl;
    //int nthreads, tid;
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
    //thrust::counting_iterator<std::size_t> ex0(0);  
    //thrust::counting_iterator<std::size_t> ex1(nthreads-1);
    //              thrust::for_each(thrust::device, ex0,ex1,
    //                               ompPrint());
  #endif

  float dt;
  const int nP = cfg.lookup("impurityParticleSource.nP");
  long nParticles = nP;
  int nT;

  if (cfg.lookupValue("timeStep.dt",dt) &&
      cfg.lookupValue("timeStep.nT",nT))    
  {
    cout << "Number of time steps: " << nT << " With dt = " << dt << endl; 
    cout << "Number of particles: " << nP << endl;              
  }
  else
  {std::cout << "ERROR: could not get nT, dt, or nP from input file" << std::endl;}

  auto particleArray = new Particles(nParticles);
  auto particleArray2 = new Particles(nParticles);
  
  float x,y,z,E,vx,vy,vz,Ex,Ey,Ez,amu,Z,charge,phi,theta,Ex_prime,Ez_prime,theta_transform;      
  if(cfg.lookupValue("impurityParticleSource.initialConditions.impurity_amu",amu) && 
     cfg.lookupValue("impurityParticleSource.initialConditions.impurity_Z",Z) &&
     cfg.lookupValue("impurityParticleSource.initialConditions.charge",charge))
    { std::cout << "Impurity amu Z charge: " << amu << " " << Z << " " << charge << std::endl;
    }
    else
    { std::cout << "ERROR: Could not get point source impurity initial conditions" << std::endl;}
  int nSourceSurfaces = 0; 
  #if PARTICLE_SOURCE_SPACE == 0 // Point Source
    if (cfg.lookupValue("impurityParticleSource.initialConditions.x_start",x) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.y_start",y) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.z_start",z))
    { std::cout << "Impurity point source: " << x << " " << y << " " << z << std::endl;
    }
    else
    { std::cout << "ERROR: Could not get point source impurity initial conditions" << std::endl;}
  #elif PARTICLE_SOURCE_SPACE > 0 //Material Surfaces - flux weighted source
    Config cfg_particles;
    std::string particleSourceFile; 
    getVariable(cfg,"particleSource.fileString",particleSourceFile);
    std::cout << "Open particle source file " << input_path+particleSourceFile << std::endl; 
    importLibConfig(cfg_particles,input_path+particleSourceFile);
    std::cout << "Successfully staged input and particle source file " << std::endl;
    
    Setting& particleSourceSetting = cfg_particles.lookup("particleSource");
    std::cout << "Successfully set particleSource setting " << std::endl;
    int nSourceBoundaries = 0,nSourceElements=0;
    float sourceMaterialZ = 0.0,accumulatedLengthArea = 0.0,sourceSampleResolution = 0.0;
    if (cfg_particles.lookupValue("particleSource.materialZ",sourceMaterialZ))
    {std::cout << "Particle Source Material Z: " << sourceMaterialZ << std::endl;}
    else
    { std::cout << "ERROR: Could not get particle source material Z" << std::endl;}
    if(sourceMaterialZ > 0.0)
    {   //count source boundaries 
        for(int i=0;i<nLines;i++)
        {
            if(boundaries[i].Z == sourceMaterialZ)
            {
                nSourceBoundaries++;
                #if USE3DTETGEOM
                  accumulatedLengthArea = accumulatedLengthArea+boundaries[i].area;
                #else
                  accumulatedLengthArea = accumulatedLengthArea+boundaries[i].length;
                #endif
            }
        }
    }
    else
    {
      if (cfg_particles.lookupValue("particleSource.nSourceBoundaries",nSourceBoundaries))
      {std::cout << "Particle Source nSourceBoundaries: " << nSourceBoundaries << std::endl;}    
      else
      { std::cout << "ERROR: Could not get particle source nSourceBoundaries" << std::endl;}
      for(int i=0;i<nSourceBoundaries;i++)
      {
          #if USE3DTETGEOM
            accumulatedLengthArea = accumulatedLengthArea+boundaries[int(particleSourceSetting["surfaceIndices"][i])].area;
          #else
            accumulatedLengthArea = accumulatedLengthArea+boundaries[int(particleSourceSetting["surfaceIndices"][i])].length;
          #endif
      }
    }
    if (cfg_particles.lookupValue("particleSource.sourceSampleResolution",sourceSampleResolution))
    {std::cout << "Particle Source sample resolution: " << sourceSampleResolution << std::endl;}
    else
    { std::cout << "ERROR: Could not get particle source sample resolution" << std::endl;}
    nSourceElements = ceil(accumulatedLengthArea/sourceSampleResolution);
    std::cout << "nSourceBoundaries accumulatedLength nSourceElements " << nSourceBoundaries << " " 
    << accumulatedLengthArea << " " << nSourceElements << std::endl;
    sim::Array<float> particleSourceSpaceCDF(nSourceElements,0.0),particleSourceX(nSourceElements,0.0),
      particleSourceY(nSourceElements,0.0),particleSourceZ(nSourceElements,0.0),
      particleSourceSpaceGrid(nSourceElements,0.0);
    sim::Array<int> particleSourceIndices(nSourceElements,0),
      particleSourceBoundaryIndices(nSourceBoundaries,0);
    #if PARTICLE_SOURCE_SPACE == 1
      for(int i=0;i<nSourceBoundaries;i++)
      {
          particleSourceBoundaryIndices[i] = particleSourceSetting["surfaceIndices"][i];
      }
      int currentSegmentIndex=0,currentBoundaryIndex=0;
      float currentAccumulatedLengthArea=0.0,lengthAlongBoundary=0.0,bDotSurfaceNorm=0.0;
      float parVec[3] = {0.0};
      float perpVec[3] = {0.0};
      currentBoundaryIndex=particleSourceBoundaryIndices[currentSegmentIndex];
      currentAccumulatedLengthArea=currentAccumulatedLengthArea+boundaries[currentBoundaryIndex].length;
      for(int i=0;i<nSourceElements;i++)
      { 
        if(i*sourceSampleResolution > currentAccumulatedLengthArea)
        {
          currentSegmentIndex++;  
          currentBoundaryIndex = particleSourceBoundaryIndices[currentSegmentIndex];
          currentAccumulatedLengthArea=currentAccumulatedLengthArea+boundaries[currentBoundaryIndex].length;  
        }
        particleSourceIndices[i] = currentBoundaryIndex;
        particleSourceBoundaryIndices[currentSegmentIndex] = particleSourceSetting["surfaceIndices"][currentSegmentIndex];
        boundaries[currentBoundaryIndex].getSurfaceParallel(parVec);
        lengthAlongBoundary = i*sourceSampleResolution - (currentAccumulatedLengthArea-boundaries[currentBoundaryIndex].length);
        particleSourceX[i]=boundaries[currentBoundaryIndex].x1 + parVec[0]*lengthAlongBoundary;
        particleSourceZ[i]=boundaries[currentBoundaryIndex].z1 + parVec[2]*lengthAlongBoundary;
        float localN = interp2dCombined(particleSourceX[i],0.0,particleSourceZ[i],nR_Dens,
                      nZ_Dens,DensGridr.data(),DensGridz.data(),ni.data());
        float localT = interp2dCombined(particleSourceX[i],0.0,particleSourceZ[i],nR_Temp,
                      nZ_Temp,TempGridr.data(),TempGridz.data(),ti.data());
        float localCs = sqrt(2*localT*1.602e-19/(1.66e-27*background_amu));
        float localBnorm[3] = {0.0}; 
          interp2dVector(&localBnorm[0],particleSourceX[i],0.0,particleSourceZ[i],nR_Bfield,
                      nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
        vectorNormalize(localBnorm,localBnorm);
        boundaries[currentBoundaryIndex].getSurfaceNormal(perpVec);
        bDotSurfaceNorm = abs(vectorDotProduct(localBnorm,perpVec));
        float localY = interp2dCombined(log10(3.0*localT),0.0,acos(bDotSurfaceNorm)*180/3.14159,nE_surfaceModel,
                      nA_surfaceModel,Elog_surfaceModel.data(),A_surfaceModel.data(),spyl_surfaceModel.data());
        localY = interp2dCombined(acos(bDotSurfaceNorm)*180/3.14159,0.0,log10(3.0*localT),nA_surfaceModel,
                      nE_surfaceModel,A_surfaceModel.data(),Elog_surfaceModel.data(),spyl_surfaceModel.data());
        std::cout << "LocalPotential localAngle localY " << 3.0*localT << " " <<acos(bDotSurfaceNorm)*180/3.1415 << " " << localY << std::endl;
        float localFlux=localCs*localN*bDotSurfaceNorm;//dotB*surf
        std::cout << "segment boundary pos x z n t cs flux " << i << " " << currentBoundaryIndex
            << " " <<particleSourceX[i] << " " << particleSourceZ[i] << " " << localN << " " << 
            localT << " " << localCs << " " <<localFlux<<std::endl;
        std::cout << "bfield perpvec bDotSurf " << localBnorm[0] << " " << localBnorm[1]
            << " " << localBnorm[2] << " " << perpVec[0] << " " << perpVec[1] << " " << 
            perpVec[2] << " " << bDotSurfaceNorm << " " << acos(bDotSurfaceNorm)*180/3.1415<<
           " "<< localY <<std::endl;
        if(i==0)
        {
          particleSourceSpaceCDF[i] = localFlux*localY; 
        }
        else
        {
          particleSourceSpaceCDF[i] = particleSourceSpaceCDF[i-1]+localFlux*localY; 
        }
          std::cout << "particleSourceSpaceCDF " << i << " " << particleSourceSpaceCDF[i] << std::endl;
      }
      for(int i=0;i<nSourceElements;i++)
      {
          particleSourceSpaceCDF[i] = particleSourceSpaceCDF[i]/particleSourceSpaceCDF[nSourceElements-1];
          std::cout << "particleSourceSpaceCDF " << i << " " << particleSourceIndices[i] << " " << 
             particleSourceX[i] << " " << particleSourceZ[i] << particleSourceSpaceCDF[i] << std::endl;
      }
      std::random_device randDevice;
      boost::random::mt19937 s0;
      s0.seed(123456);
      boost::random::uniform_01<> dist01;
      float rand0 = 0.0;
      int lowInd = 0;
      int currentSegment = 0;
    #else
    #endif
  #endif    
  #if PARTICLE_SOURCE_ENERGY == 0
    if( cfg.lookupValue("impurityParticleSource.initialConditions.energy_eV",E))
    { std::cout << "Impurity point source E: " << E << std::endl;
    }
    else
    { std::cout << "ERROR: Could not get point source impurity initial conditions" << std::endl;}
  #elif PARTICLE_SOURCE_ENERGY > 0
    #if PARTICLE_SOURCE_ENERGY == 1
    //Create Thompson Distribution
      float surfaceBindingEnergy = cfg.lookup("impurityParticleSource.source_material_SurfaceBindingEnergy");
      float surfaceAlpha = cfg.lookup("impurityParticleSource.source_materialAlpha");
      std::cout << "surface binding energy " << surfaceBindingEnergy << std::endl;
      int nThompDistPoints = 200;
      float max_Energy = 100.0;
      sim::Array<float> ThompsonDist(nThompDistPoints),CumulativeDFThompson(nThompDistPoints);
      for(int i=0;i<nThompDistPoints;i++)
      {
        if(surfaceAlpha > 0.0)
        {
          ThompsonDist[i] = surfaceAlpha*(surfaceAlpha-1.0)*(i*max_Energy/nThompDistPoints)*pow(surfaceBindingEnergy,surfaceAlpha-1.0)/pow((i*max_Energy/nThompDistPoints) + surfaceBindingEnergy,(surfaceAlpha+1.0));
        }
        else
        {  
          ThompsonDist[i] = (i*max_Energy/nThompDistPoints)/pow((i*max_Energy/nThompDistPoints) + surfaceBindingEnergy,3);
        }
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
      #elif PARTICLE_SOURCE_ENERGY == 2
      #endif
        boost::random::mt19937 sE;
        boost::random::uniform_01<> dist01E;
        float randE = 0.0;
        int lowIndE = 0;
  #endif
  #if PARTICLE_SOURCE_ANGLE == 0
    if (cfg.lookupValue("impurityParticleSource.initialConditions.phi",phi) &&
        cfg.lookupValue("impurityParticleSource.initialConditions.theta",theta))
    { std::cout << "Impurity point source angles phi theta: " << phi << " " << theta << std::endl;
    }
    else
    { std::cout << "ERROR: Could not get point source angular initial conditions" << std::endl;}
    phi = phi*3.141592653589793/180.0;
    theta = theta*3.141592653589793/180.0;
    Ex = E*sin(phi)*cos(theta);
    Ey = E*sin(phi)*sin(theta);
    Ez = E*cos(phi);
    if(phi == 0.0)
    {
        Ex=0.0;
        Ey=0.0;
        Ez=E;
    }
  #elif PARTICLE_SOURCE_ANGLE > 0

    std::cout << "Read particle source " << std::endl;
    #if PARTICLE_SOURCE_ENERGY < 2
      Config cfg_particles;
    #endif
    //cfg_particles.readFile((input_path+"particleSource.cfg").c_str());
    //Setting& particleSource = cfg_particles.lookup("particleSource");
    //int nSegmentsAngle = particleSource["nSegmentsAngle"];
    //float angleSample;
    //sim::Array<float> sourceAngleSegments(nSegmentsAngle);
    //sim::Array<float> angleCDF(nSegmentsAngle);
    //for (int i=0; i<(nSegmentsAngle); i++)
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
  #if PARTICLE_SOURCE_FILE > 0 // File source
    Config cfg_particles;
    std::string ncParticleSourceFile;
    getVariable(cfg,"particleSource.ncFileString",ncParticleSourceFile);
std::cout << "About to try to open NcFile ncp0 " << std::endl;
// Return this in event of a problem.
static const int NC_ERR = 2;   
try
   {
    NcFile ncp0("input/"+ncParticleSourceFile, NcFile::read);
   }catch(NcException& e)
     {
       e.what();
       cout<<"FAILURE*************************************"<<endl;
       return NC_ERR;
     }
std::cout << "finished NcFile ncp0 starting ncp" << std::endl;
    NcFile ncp("input/"+ncParticleSourceFile, NcFile::read);
std::cout << "getting dim nP" << std::endl;
    NcDim ps_nP(ncp.getDim("nP"));

    int nPfile = ps_nP.getSize();

std::cout << "nPfile "<< nPfile << std::endl;
    NcVar ncp_x(ncp.getVar("x"));
    NcVar ncp_y(ncp.getVar("y"));
    NcVar ncp_z(ncp.getVar("z"));
    NcVar ncp_vx(ncp.getVar("vx"));
    NcVar ncp_vy(ncp.getVar("vy"));
    NcVar ncp_vz(ncp.getVar("vz"));

std::cout << "got through NcVar " << std::endl;
    vector<float> xpfile(nPfile),ypfile(nPfile),zpfile(nPfile),
    vxpfile(nPfile),vypfile(nPfile),vzpfile(nPfile);
std::cout << "defined file vectors " << std::endl;
    ncp_x.getVar(&xpfile[0]);  
    ncp_y.getVar(&ypfile[0]);  
    ncp_z.getVar(&zpfile[0]);  
    ncp_vx.getVar(&vxpfile[0]);  
    ncp_vy.getVar(&vypfile[0]);  
    ncp_vz.getVar(&vzpfile[0]);
std::cout << "defined file vectors " << std::endl;
  ncp.close();  
std::cout << "closed ncp " << std::endl;
    //for(int i=0;i<nPfile;i++)
    //{
    //    std::cout << " xyz from file " << xpfile[i] << " " << ypfile[i] << " " << zpfile[i] << std::endl;  
    //    std::cout << " Exyz from file " << Expfile[i] << " " << Eypfile[i] << " " << Ezpfile[i] << std::endl;  
    //}
  #endif
  sim::Array<float> pSurfNormX(nP),pSurfNormY(nP),pSurfNormZ(nP), 
                    px(nP),py(nP),pz(nP),pvx(nP),pvy(nP),pvz(nP);
  int surfIndexMod = 0;
  float eVec[3] = {0.0};
  for (int i=0; i< nP ; i++)
  {
    #if PARTICLE_SOURCE_SPACE > 0 // File source
      #if USE3DTETGEOM > 0
      surfIndexMod = i%nSourceSurfaces;
      float xCentroid = (boundaries[sourceElements[surfIndexMod]].x1 + boundaries[sourceElements[surfIndexMod]].x2 + boundaries[sourceElements[surfIndexMod]].x3)/3.0;
      float yCentroid = (boundaries[sourceElements[surfIndexMod]].y1 + boundaries[sourceElements[surfIndexMod]].y2 + boundaries[sourceElements[surfIndexMod]].y3)/3.0;
      float zCentroid = (boundaries[sourceElements[surfIndexMod]].z1 + boundaries[sourceElements[surfIndexMod]].z2 + boundaries[sourceElements[surfIndexMod]].z3)/3.0;
      float bufferLaunch = 1.0e-4;
      x = xCentroid - bufferLaunch*boundaries[sourceElements[surfIndexMod]].a/boundaries[sourceElements[surfIndexMod]].plane_norm;//boundaries[sourceElements[surfIndexMod]].x1;
      y = yCentroid - bufferLaunch*boundaries[sourceElements[surfIndexMod]].b/boundaries[sourceElements[surfIndexMod]].plane_norm;//boundaries[sourceElements[surfIndexMod]].y1;
      z = zCentroid - bufferLaunch*boundaries[sourceElements[surfIndexMod]].c/boundaries[sourceElements[surfIndexMod]].plane_norm;//boundaries[sourceElements[surfIndexMod]].z1; 
      #else
        //x = sampled
        rand0 = dist01(s0);
        float distAlongSegs = interp1dUnstructured(rand0,nSourceElements,accumulatedLengthArea,&particleSourceSpaceCDF[0], lowInd);
        currentSegment = particleSourceIndices[lowInd];
        std::cout << "rand of " << rand0 << " puts the particle " << distAlongSegs << " along the segments on the boundary element " << currentSegment << std::endl;
        float parVec[3] = {0.0};
        boundaries[currentSegment].getSurfaceParallel(parVec);
        x = particleSourceX[lowInd]+(rand0-particleSourceSpaceCDF[lowInd])/(particleSourceSpaceCDF[lowInd+1]-particleSourceSpaceCDF[lowInd])*sourceSampleResolution*parVec[0];
        y = 0.0;
        z = particleSourceZ[lowInd]+(rand0-particleSourceSpaceCDF[lowInd])/(particleSourceSpaceCDF[lowInd+1]-particleSourceSpaceCDF[lowInd])*sourceSampleResolution*parVec[2];
        float buffer = 1e-6;//0.0;//2e-6;
        x = x - buffer*boundaries[currentSegment].a/boundaries[currentSegment].plane_norm;//boundaries[sourceElements[surfIndexMod]].x1;
        z = z - buffer*boundaries[currentSegment].c/boundaries[currentSegment].plane_norm;//boundaries[sourceElements[surfIndexMod]].z1; 
      #endif
    #endif
    #if PARTICLE_SOURCE_ENERGY > 0
        randE = dist01E(sE);
      #if PARTICLE_SOURCE_ENERGY == 1
        E = interp1dUnstructured(randE,nThompDistPoints, max_Energy, &CumulativeDFThompson.front(), lowIndE);
      #elif PARTICLE_SOURCE_ENERGY == 2
        float localT = interp2dCombined(x,y,z,nR_Temp,
                      nZ_Temp,TempGridr.data(),TempGridz.data(),ti.data());
        float localBnorm[3] = {0.0}; 
          interp2dVector(&localBnorm[0],x,y,z,nR_Bfield,
                      nZ_Bfield,bfieldGridr.data(),bfieldGridz.data(),br.data(),bz.data(),by.data());
        vectorNormalize(localBnorm,localBnorm);
        boundaries[currentSegment].getSurfaceNormal(perpVec);
        bDotSurfaceNorm = abs(vectorDotProduct(localBnorm,perpVec));
        float localAngle = acos(bDotSurfaceNorm)*180/3.1415;
        float sputtE = interp3d ( randE,localAngle,log10(3.0*localT),nEdistBins_surfaceModel,nA_surfaceModel,nE_surfaceModel,
              energyDistGrid01.data(),A_surfaceModel.data(),Elog_surfaceModel.data() ,energyDist_CDFregrid.data() );
        E = sputtE;
        std::cout << "randE of " << randE << " with localAngle " << localAngle << " and local potential " <<
           3.0*localT << " puts the particle energy to " << E  << std::endl;
      #endif  
    #endif    
    #if PARTICLE_SOURCE_ANGLE == 1 // Analytic normal incidence
      Ex = -E*boundaries[currentSegment].a/boundaries[currentSegment].plane_norm;
      Ey = -E*boundaries[currentSegment].b/boundaries[currentSegment].plane_norm;
      Ez = -E*boundaries[currentSegment].c/boundaries[currentSegment].plane_norm;
    
    #elif PARTICLE_SOURCE_ANGLE > 1
      randA = dist01A(sA);
      float sputtA = interp3d ( randA,localAngle,log10(3.0*localT),nAdistBins_surfaceModel,nA_surfaceModel,nE_surfaceModel,
              cosDistGrid01.data(),A_surfaceModel.data(),Elog_surfaceModel.data() ,cosDist_CDFregrid.data() );
      phi = sputtA*3.141592653589793/180.0;
      std::cout << "sputtA and phi " << sputtA << " " << phi << std::endl;
      randA = dist01A(sA);
      theta = 2.0*3.141592653589793*randA;
      std::cout << "randA and theta " << randA << " " << theta << std::endl;
      Ex = E*sin(phi)*cos(theta);
      Ey = E*sin(phi)*sin(theta);
      Ez = E*cos(phi);
      std::cout << "randA of " << randA << " puts the particle angle phi to " << phi  << std::endl;
      std::cout << "E of particle " << Ex << " " << Ey << " " << Ez << " " << std::endl;
      std::cout << "current segment and perpVec " << currentSegment << " " << perpVec[0] << " " << perpVec[1] << " " << perpVec[2] << std::endl;
      float Ezx = sqrt(Ez*Ez + Ex*Ex);
      float thetaEzx = atan2(Ez,Ex);
      std::cout << "Ezx thetaEzx " << Ezx << " " << thetaEzx << std::endl;
      //positive slope equals negative upward normal
      theta_transform = acos(perpVec[2]);//-sgn(boundaries[currentSegment].slope_dzdx)*
      //if(perpVec[2]==0.0)
      //{
      //    if(perpVec[0] > 0.0)
      //    {
      //      theta_transform = 0.5*3.141592653589793;
      //      std::cout << "Vertical line element perpVec " << perpVec[0] << " " << perpVec[1] << " " << perpVec[2] << " " << theta_transform << std::endl;
      //    }
      //    else if(perpVec[0] < 0.0)
      //    {
      //      theta_transform = 1.5*3.141592653589793;
      //      std::cout << "Vertical line element perpVec " << perpVec[0] << " " << perpVec[1] << " " << perpVec[2] << " " << theta_transform << std::endl;
      //    }
      //}
      Ex = Ezx*cos(thetaEzx - theta_transform);
      //Ey = E*sin(phi+theta_transform)*sin(theta);
      Ez = Ezx*sin(thetaEzx - theta_transform);
      std::cout << "theta transform " << theta_transform << std::endl;
      eVec[0] = Ex;
      eVec[1] = Ey;
      eVec[2] = Ez;
      float EdotP = vectorDotProduct(perpVec,eVec);
      if(EdotP < 0.0)
      {
          std::cout << "This dot product negative " << std::endl;
          Ex = -Ex;
          Ez = -Ez;
      }
      //Ex_prime = Ex*cos(theta_transform) - Ez*sin(theta_transform);
      //Ez_prime = Ex*sin(theta_transform) + Ez*cos(theta_transform);
      //Ex = Ex_prime;
      //Ez = Ez_prime;
      std::cout << "Transformed E " << Ex << " " << Ey << " " << Ez << " " << std::endl;
      //particleArray->setParticle(i,x, y, z, Ex, Ey,Ez, Z, amu, charge);
    #endif

    #if PARTICLE_SOURCE_FILE > 0 // File source
      x = xpfile[i];
      y = ypfile[i];
      z = zpfile[i];
      vx = vxpfile[i];
      vy = vypfile[i];
      vz = vzpfile[i];
    #endif  
//    std::cout << "particle xyz Exyz Z amu charge " << x << " " << y << " " << z << " "
//       << Ex << " " << Ey << " " << Ez << " " << Z << " " << amu << " " << charge << " "  << std::endl;
    particleArray->setParticleV(i,x,y, z, vx, vy, vz, Z, amu, charge);   
    #if PARTICLE_SOURCE_SPACE > 0
    pSurfNormX[i] = -boundaries[currentSegment].a/boundaries[currentSegment].plane_norm;
    pSurfNormY[i] = -boundaries[currentSegment].b/boundaries[currentSegment].plane_norm;
    pSurfNormZ[i] = -boundaries[currentSegment].c/boundaries[currentSegment].plane_norm;
    #endif
    px[i] = x;
    py[i] = y;
    pz[i] = z;
    pvx[i] = vx;
    pvy[i] = vy;
    pvz[i] = vz;
  } 
#if USE_MPI > 0
  if(world_rank==0)
  { 
#endif
std::cout <<" about to write ncFile_particles " << std::endl;
    NcFile ncFile_particles("output/particleSource.nc", NcFile::replace);
    std::cout <<" opened file " << std::endl;
    NcDim pNP = ncFile_particles.addDim("nP",nP);
    NcVar p_surfNormx = ncFile_particles.addVar("surfNormX",ncFloat,pNP);
    NcVar p_surfNormy = ncFile_particles.addVar("surfNormY",ncFloat,pNP);
    NcVar p_surfNormz = ncFile_particles.addVar("surfNormZ",ncFloat,pNP);
    NcVar p_vx = ncFile_particles.addVar("vx",ncFloat,pNP);
    NcVar p_vy = ncFile_particles.addVar("vy",ncFloat,pNP);
    NcVar p_vz = ncFile_particles.addVar("vz",ncFloat,pNP);
    NcVar p_x = ncFile_particles.addVar("x",ncFloat,pNP);
    NcVar p_y = ncFile_particles.addVar("y",ncFloat,pNP);
    NcVar p_z = ncFile_particles.addVar("z",ncFloat,pNP);
    std::cout <<" added vars " << std::endl;
    p_surfNormx.putVar(&pSurfNormX[0]);
    p_surfNormy.putVar(&pSurfNormY[0]);
    p_surfNormz.putVar(&pSurfNormZ[0]);
    p_vx.putVar(&pvx[0]);
    p_vy.putVar(&pvy[0]);
    p_vz.putVar(&pvz[0]);
    p_x.putVar(&px[0]);
    p_y.putVar(&py[0]);
    p_z.putVar(&pz[0]);
    std::cout <<" put vars complete " << std::endl;
    ncFile_particles.close();
NcFile ncFile_test("testfile.nc", NcFile::replace);
ncFile_test.close();
std::cout <<" closed ncFile_particles " << std::endl;
#if USE_MPI > 0
  }
#endif

  std::cout << "finished loading particle source" << std::endl;

  #if GEOM_TRACE > 0       
    std::uniform_real_distribution<float> dist2(0,1);
    //std::random_device rd2;
    //std::default_random_engine generator2(rd2());
    float randDevice02 = 6.52E+5;
    std::default_random_engine generatorTrace(randDevice02);
    std::cout << "Randomizing velocities to trace geometry. " << std::endl;

    for (int i=0 ; i<nParticles ; i++)
    {   float theta_trace = dist2(generatorTrace)*2*3.1415;
        float phi_trace = dist2(generatorTrace)*3.1415;
        float mag_trace = 2e3;
        particleArray->vx[i] = mag_trace*cos(theta_trace)*sin(phi_trace);
        particleArray->vy[i] = mag_trace*sin(theta_trace)*sin(phi_trace);
        particleArray->vz[i] = mag_trace*cos(phi_trace);
    }
  #endif

  #if PARTICLE_TRACKS > 0
    int subSampleFac = 1;
    if(cfg.lookupValue("diagnostics.trackSubSampleFactor", subSampleFac))
       {std::cout << "Tracks subsample factor imported" << std::endl;}
    else
    { std::cout << "ERROR: Could not get tracks sub sample info from input file " << std::endl;}
    std::cout << "history array length " << (nT/subSampleFac)*nP << std::endl;
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
      positionHistoryX[0] = new float [(nT/subSampleFac)*nP];
      positionHistoryY[0] = new float [(nT/subSampleFac)*nP];
      positionHistoryZ[0] = new float [(nT/subSampleFac)*nP];
      velocityHistoryX[0] = new float [(nT/subSampleFac)*nP];
      velocityHistoryY[0] = new float [(nT/subSampleFac)*nP];
      velocityHistoryZ[0] = new float [(nT/subSampleFac)*nP];
      chargeHistory[0] = new float [(nT/subSampleFac)*nP];
      for(int i=0 ; i<nP ; i++)
      {
          positionHistoryX[i] = &positionHistoryX[0][(nT/subSampleFac)*i];
          positionHistoryY[i] = &positionHistoryY[0][(nT/subSampleFac)*i];
          positionHistoryZ[i] = &positionHistoryZ[0][(nT/subSampleFac)*i];
          velocityHistoryX[i] = &velocityHistoryX[0][(nT/subSampleFac)*i];
          velocityHistoryY[i] = &velocityHistoryY[0][(nT/subSampleFac)*i];
          velocityHistoryZ[i] = &velocityHistoryZ[0][(nT/subSampleFac)*i];
          chargeHistory[i] = &chargeHistory[0][(nT/subSampleFac)*i];
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

  std::cout << "beginning seeds" << std::endl;
  std::uniform_real_distribution<float> dist(0,1e6);

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

  thrust::counting_iterator<std::size_t> particleBegin(0);  
  thrust::counting_iterator<std::size_t> particleEnd(nParticles);
    
  #if PARTICLESEEDS > 0
    #if USE_CUDA  
      sim::Array<curandState> state1(nParticles);
    #else
      sim::Array<std::mt19937> state1(nParticles);
    #endif
    #if USEIONIZATION > 0 || USERECOMBINATION > 0 || USEPERPDIFFUSION > 0 || USECOULOMBCOLLISIONS > 0 || USESURFACEMODEL > 0
      std::cout << "Initializing curand seeds " << std::endl;
      thrust::for_each(thrust::device, particleBegin,particleEnd,
                           curandInitialize(&state1[0],0));
      #if USE_CUDA
        cudaDeviceSynchronize();
      #endif
    #endif
  #endif

    float moveTime = 0.0;
    float geomCheckTime = 0.0;
    float ionizTime = 0.0;

#if USE_BOOST
    //cpu_times copyToDeviceTime = timer.elapsed();
    //std::cout << "Initialize rand state and copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';
#endif
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::duration<float> fsec;
    auto start_clock = Time::now();
    std::cout << "Starting main loop"  << std::endl;
    float testFlowVec[3] = {0.0f};
    interp2dVector(&testFlowVec[0],0.01,-0.02,0.1,nR_flowV,nZ_flowV,
                                     flowVGridr.data(),flowVGridz.data(),
                                     flowVr.data(),flowVz.data(),flowVt.data());    
std::cout << "Flow vNs "<< testFlowVec[0] << " " <<testFlowVec[1] << " " << testFlowVec[2]  << std::endl;
    std::cout << "Starting main loop" << particleArray->xprevious[0] << std::endl;
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
#endif
    std::cout << "nDevices " << nDevices  << "  " << cpu_thread_id << " " << num_cpu_threads<< " particle index " << cpu_thread_id*nP/nDevices << " " << (cpu_thread_id+1)*nP/nDevices - 1 << std::endl;
//int cpu_thread_id = 0;
//nDevices = 2;
#endif
    for(int tt=0; tt< nT; tt++)
    {
#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif
        thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleEnd, 
                move_boris(particleArray,dt,boundaries.data(), nLines,
                    nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&by.front(),
                    nR_PreSheathEfield,nY_PreSheathEfield,nZ_PreSheathEfield,
                    &preSheathEGridr.front(),&preSheathEGridy.front(),&preSheathEGridz.front(),
                    &PSEr.front(),&PSEz.front(),&PSEt.front(),
                        nR_closeGeom_sheath,nY_closeGeom_sheath,nZ_closeGeom_sheath,n_closeGeomElements_sheath,
                        &closeGeomGridr_sheath.front(),&closeGeomGridy_sheath.front(),&closeGeomGridz_sheath.front(),
                        &closeGeom_sheath.front()) );
        //std::cout << "pos " << particleArray->x[0] << "  " << particleArray->xprevious[0] << std::endl;   
       //if(std::isnan(particleArray->x[0]))
         // exit(0);
       //particleArray->x[0] = 0.0;
        //try {
            thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin, particleEnd,
                    geometry_check(particleArray,nLines,&boundaries[0],surfaces,dt,tt,
                        nR_closeGeom,nY_closeGeom,nZ_closeGeom,n_closeGeomElements,
                        &closeGeomGridr.front(),&closeGeomGridy.front(),&closeGeomGridz.front(),
                        &closeGeom.front()) );
       // }
       /*
            catch (thrust::system_error &e) {
            std::cerr << "Thrust system error: " << e.what() << std::endl;
            exit(-1);
        }
        */
#if SPECTROSCOPY > 0
            thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,// particleBegin,particleEnd,
                    spec_bin(particleArray,nBins,net_nX,net_nY, net_nZ, &gridX_bins.front(),&gridY_bins.front(),
                        &gridZ_bins.front(), &net_Bins.front(),dt) );
#endif            
#if USEIONIZATION > 0
        thrust::for_each(thrust::device, particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin,particleEnd,
                ionize(particleArray, dt,&state1.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),  
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                    nTemperaturesIonize, nDensitiesIonize,&gridTemperature_Ionization.front(),
                    &gridDensity_Ionization.front(), &rateCoeff_Ionization.front(),tt));
#endif
#if USERECOMBINATION > 0
        thrust::for_each(thrust::device, particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin,particleEnd,
                recombine(particleArray, dt,&state1.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),  
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                    nTemperaturesRecombine,nDensitiesRecombine,
                    gridTemperature_Recombination.data(),gridDensity_Recombination.data(),
                    rateCoeff_Recombination.data(),tt));
#endif
#if USEPERPDIFFUSION > 0
        thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin, particleEnd,
                crossFieldDiffusion(particleArray,dt,&state1.front(),perpDiffusionCoeff,
                    nR_Bfield,nZ_Bfield,bfieldGridr.data(),&bfieldGridz.front(),
                                        &br.front(),&bz.front(),&by.front()));
            
            thrust::for_each(thrust::device, particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin,particleEnd,
                    geometry_check(particleArray,nLines,&boundaries[0],surfaces,dt,tt,
                        nR_closeGeom,nY_closeGeom,nZ_closeGeom,n_closeGeomElements,
                        &closeGeomGridr.front(),&closeGeomGridy.front(),&closeGeomGridz.front(),
                        &closeGeom.front()) );
#endif
#if USECOULOMBCOLLISIONS > 0
        thrust::for_each(thrust::device, particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin, particleEnd, 
                coulombCollisions(particleArray,dt,&state1.front(),
                    nR_flowV,nY_flowV,nZ_flowV,&flowVGridr.front(),&flowVGridy.front(),&flowVGridz.front(),
                    &flowVr.front(),&flowVz.front(),&flowVt.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),    
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),ti.data(),&te.front(),
                    background_Z,background_amu, 
                    nR_Bfield,nZ_Bfield,bfieldGridr.data(),&bfieldGridz.front(),
                                        &br.front(),&bz.front(),&by.front()));

#endif
#if USETHERMALFORCE > 0
        thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin, particleEnd,
                thermalForce(particleArray,dt,background_amu,
                    nR_gradT,nZ_gradT,gradTGridr.data(),gradTGridz.data(),
                    gradTiR.data(),gradTiZ.data(), gradTiY.data(), 
                    gradTeR.data(), gradTeZ.data(), gradTeY.data(), 
                    nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&by.front()));
#endif

#if USESURFACEMODEL > 0
        thrust::for_each(thrust::device,particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin, particleEnd, 
                reflection(particleArray,dt,&state1.front(),nLines,&boundaries[0],surfaces,
                    nE_sputtRefCoeff, nA_sputtRefCoeff,A_sputtRefCoeff.data(),
                    Elog_sputtRefCoeff.data(),spyl_surfaceModel.data(), rfyl_surfaceModel.data(),
                    nE_sputtRefDistOut, nA_sputtRefDistOut,nE_sputtRefDistIn,nA_sputtRefDistIn,
                    Elog_sputtRefDistIn.data(),A_sputtRefDistIn.data(),
                    E_sputtRefDistOut.data(),A_sputtRefDistOut.data(),
                    energyDistGrid01.data(),angleDistGrid01.data(),
                    EDist_CDF_Y_regrid.data(),ADist_CDF_Y_regrid.data(),
                    EDist_CDF_R_regrid.data(),ADist_CDF_R_regrid.data(),
                    nEdist, E0dist, Edist, nAdist, A0dist, Adist) );
#endif        

#if PARTICLE_TRACKS >0
#if USE_CUDA > 0
   thrust::for_each(thrust::device, particleBegin+ world_rank*nP/world_size,particleBegin + (world_rank+1)*nP/world_size,//particleBegin,particleEnd,
      history(particleArray,tt,subSampleFac,nP,&positionHistoryX.front(),
      &positionHistoryY.front(),&positionHistoryZ.front(),
      &velocityHistoryX.front(),&velocityHistoryY.front(),
      &velocityHistoryZ.front(),&chargeHistory.front()) );
#else
if (tt % subSampleFac == 0)  
{    
        for(int i=0;i<nP;i++)
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
    }
#if USE_OPENMP
}
#endif
// Ensure that all time step loop GPU kernels are complete before proceeding
    #ifdef __CUDACC__
        cudaDeviceSynchronize();
    #endif

    auto finish_clock = Time::now();
    fsec fs = finish_clock - start_clock;
    printf("Time taken          is %6.3f (secs) \n", fs.count());
    printf("Time taken per step is %6.3f (secs) \n", fs.count() / (float) nT);
#if USE_BOOST
    //cpu_times ionizeTimeGPU = timer.elapsed();
    //std::cout << "Particle Moving Time: " << ionizeTimeGPU.wall*1e-9 << '\n';
#endif
//for(int i=0; i<nP;i++)
//{
//    std::cout << "Particle test value r1: " << i << " " << particleArray->test[i] << std::endl;
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
    std::cout << "transit time counting "<< nP << " " << particleArray->x[0] <<  std::endl;
    //float tmp202 =0.0;
#if USE_CUDA
    cudaDeviceSynchronize();
#endif
#if USE_MPI > 0
    sim::Array<float> xGather(nP,0);
    //float *x_gather = NULL;
    //if (world_rank == 0) {
    //      x_gather = malloc(sizeof(float)*nP);
    //}
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&particleArray->x[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, &xGather[0], nP/world_size,MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    //Collect stuff
   for(int rr=1; rr<world_size;rr++)
{
if(world_rank == rr)
{
    //MPI_Send(&particleArray->x[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->y[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->z[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->vx[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->vy[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->vz[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->hitWall[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->weight[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&particleArray->charge[world_rank*nP/world_size], nP/world_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
}
else if(world_rank == 0)
{
    //MPI_Recv(&particleArray->x[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->y[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->z[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->vx[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->vy[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->vz[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->hitWall[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->weight[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particleArray->charge[rr*nP/world_size], nP/world_size, MPI_FLOAT, rr, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
#if SPECTROSCOPY > 0
//MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&net_Bins[0], &net_BinsTotal[0], (nBins+1)*net_nX*net_nZ, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
#endif
#if USESURFACEMODEL > 0
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->grossDeposition[0], &grossDeposition[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->grossErosion[0], &grossErosion[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sumWeightStrike[0], &sumWeightStrike[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->aveSputtYld[0], &aveSputtYld[0],nSurfaces, MPI_FLOAT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sputtYldCount[0], &sputtYldCount[0],nSurfaces, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->sumParticlesStrike[0], &sumParticlesStrike[0],nSurfaces, MPI_INT, MPI_SUM, 0,
                   MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
#if FLUX_EA > 0 
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(&surfaces->energyDistribution[0], &energyDistribution[0],nSurfaces*nEdist*nAdist, 
        MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
#endif
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
        std::cout << "loop " << i << std::endl;
        if(particleArray->hitWall[i] == 1.0)
        {
            meanTransitTime0 = meanTransitTime0 + particleArray->transitTime[i];
        }
    }
    */
meanTransitTime0 = meanTransitTime0/nP;
std::cout << " mean transit time " << meanTransitTime0 << std::endl;
    int max_boundary = 0;
    float max_impacts = 0.0;
    int max_boundary1 = 0;
    float max_impacts1 = 0.0;
    std::cout << " new pointers with nLines " << nLines << std::endl;
    float* impacts = new float[nLines];
    float* xOut = new float[nP];
    std::cout << " first one worked "  << std::endl;
    float* redeposit = new float[nLines];
    float* startingParticles = new float[nLines];
    float* surfZ = new float[nLines];
    //int nA = 90;
    //int nE = 1000;
    //float* impactEnergy = new float[nLines*nA*nE];
std::cout << "before starting loop "<< particleArray->xprevious[0]  << std::endl;
std::cout << " starting loop "  << std::endl;
    for (int i=0; i<nLines; i++)
    {
        impacts[i] = boundaries[i].impacts;
        redeposit[i] = boundaries[i].redeposit;
        startingParticles[i] = boundaries[i].startingParticles;
        if (boundaries[i].impacts > max_impacts)
        {
            max_impacts = boundaries[i].impacts;
            max_boundary = i;
        }
        surfZ[i] = boundaries[i].Z;
    }

    for (int i=0; i<nP; i++)
    {
xOut[i] = particleArray->x[i];
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
    float* startingParticles = new float[nLines];
    float* surfZ = new float[nLines];
    //float* impactEnergy = new float[nLines*1000];
    for (int i=0; i<nLines; i++)
    {
        impacts[i] = boundaries[i].impacts;
        startingParticles[i] = boundaries[i].startingParticles;
        surfZ[i] = boundaries[i].Z;
    }
#endif

 //add initial particle erosion to surface counting
 int closestBoundaryIndex=0;
 int surfIndex=0;
 float minDistance = 0.0;
 float thisE[3] = {0.0f};
 for(int j=0; j<nP; j++)
 {
    minDistance = getE(px[j],py[j],pz[j],thisE, 
            boundaries.data(),nLines,
            nR_closeGeom_sheath,nY_closeGeom_sheath,
            nZ_closeGeom_sheath,n_closeGeomElements_sheath,
            &closeGeomGridr_sheath.front(),&closeGeomGridy_sheath.front(),
            &closeGeomGridz_sheath.front(),&closeGeom_sheath.front(),
            closestBoundaryIndex);
    if(boundaries[closestBoundaryIndex].Z > 0.0)
    {
      surfIndex = boundaries[closestBoundaryIndex].surfaceNumber;
      grossErosion[surfIndex] = grossErosion[surfIndex] + 1.0;
    }
 }
//#if PARTICLE_SOURCE == 1
//int ring1 = 0;
//int ring2 = 0;
//int noWall = 0;
//float meanTransitTime = 0.0;
//
//for(int i=0; i<nP ; i++)
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
//meanTransitTime = meanTransitTime/(nP-noWall);
//std::cout << "Number of impurity particles deposited on ring 1 " << ring1 << std::endl;
//std::cout << "Number of impurity particles deposited on ring 2 " << ring2 << std::endl;
//std::cout << "Number of impurity particles not deposited " << noWall << std::endl;
//std::cout << "Mean transit time of deposited particles " << meanTransitTime << std::endl;
//#endif
    std::cout << "positions.m writing " << std::endl;
    ofstream outfile2;
    outfile2.open ("output/positions.m");
    for(int i=1 ; i<nP+1 ; i++)
      {
        outfile2 << "Pos( " << i<< ",:) = [ " ;
        outfile2 << particleArray->x[i-1] << " " << particleArray->y[i-1] 
            << " " << particleArray->z[i-1] << " ];" << std::endl;
      }
       outfile2.close();
       std::cout << "finished writing positions.m " << std::endl;
// Write netCDF output for positions
NcFile ncFile0("output/positions.nc", NcFile::replace);
       std::cout << "created file " << std::endl;
NcDim nc_nP0 = ncFile0.addDim("nP",nP);
       std::cout << "created dim nP " << std::endl;
vector<NcDim> dims0;
dims0.push_back(nc_nP0);
       std::cout << "created dims Vector " << std::endl;

NcVar nc_x0 = ncFile0.addVar("x",ncDouble,dims0);
NcVar nc_y0 = ncFile0.addVar("y",ncDouble,dims0);
NcVar nc_z0 = ncFile0.addVar("z",ncDouble,dims0);
NcVar nc_vx0 = ncFile0.addVar("vx",ncDouble,dims0);
NcVar nc_vy0 = ncFile0.addVar("vy",ncDouble,dims0);
NcVar nc_vz0 = ncFile0.addVar("vz",ncDouble,dims0);
NcVar nc_trans0 = ncFile0.addVar("transitTime",ncDouble,dims0);
NcVar nc_impact0 = ncFile0.addVar("hitWall",ncDouble,dims0);
NcVar nc_weight0 = ncFile0.addVar("weight",ncDouble,dims0);
NcVar nc_charge0 = ncFile0.addVar("charge",ncDouble,dims0);
       std::cout << "added Vars " << std::endl;
       std::cout << "x0 "<< particleArray->x[0] << std::endl;

nc_x0.putVar(&xOut[0]);
       std::cout << "added x " << std::endl;
nc_y0.putVar(&particleArray->y[0]);
nc_z0.putVar(&particleArray->z[0]);
nc_vx0.putVar(&particleArray->vx[0]);
nc_vy0.putVar(&particleArray->vy[0]);
nc_vz0.putVar(&particleArray->vz[0]);
nc_trans0.putVar(&particleArray->transitTime[0]);
nc_impact0.putVar(&particleArray->hitWall[0]);
nc_weight0.putVar(&particleArray->weight[0]);
nc_charge0.putVar(&particleArray->charge[0]);
ncFile0.close();
       std::cout << "closed positions opening surface " << std::endl;
#if USESURFACEMODEL > 0
NcFile ncFile1("output/surface.nc", NcFile::replace);
NcDim nc_nLines = ncFile1.addDim("nSurfaces",nSurfaces);
vector<NcDim> dims1;
dims1.push_back(nc_nLines);

vector<NcDim> dimsSurfE;
dimsSurfE.push_back(nc_nLines);
NcDim nc_nEnergies = ncFile1.addDim("nEnergies",nEdist);
NcDim nc_nAngles = ncFile1.addDim("nAngles",nAdist);
dimsSurfE.push_back(nc_nAngles);
dimsSurfE.push_back(nc_nEnergies);
NcVar nc_grossDep = ncFile1.addVar("grossDeposition",ncDouble,nc_nLines);
NcVar nc_grossEro = ncFile1.addVar("grossErosion",ncDouble,nc_nLines);
NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl",ncDouble,nc_nLines);
NcVar nc_spylCounts = ncFile1.addVar("spylCounts",ncInt,nc_nLines);
NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike",ncInt,nc_nLines);
NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike",ncDouble,nc_nLines);
nc_grossDep.putVar(&grossDeposition[0]);
nc_grossEro.putVar(&grossErosion[0]);
nc_aveSpyl.putVar(&aveSputtYld[0]);
nc_spylCounts.putVar(&sputtYldCount[0]);
nc_spylCounts.putVar(&sumParticlesStrike[0]);
nc_sumWeightStrike.putVar(&sumWeightStrike[0]);
//NcVar nc_surfImpacts = ncFile1.addVar("impacts",ncDouble,dims1);
//NcVar nc_surfRedeposit = ncFile1.addVar("redeposit",ncDouble,dims1);
//NcVar nc_surfStartingParticles = ncFile1.addVar("startingParticles",ncDouble,dims1);
//NcVar nc_surfZ = ncFile1.addVar("Z",ncDouble,dims1);
NcVar nc_surfEDist = ncFile1.addVar("surfEDist",ncDouble,dimsSurfE);
//nc_surfImpacts.putVar(impacts);
//#if USE3DTETGEOM > 0
//nc_surfRedeposit.putVar(redeposit);
//#endif
//nc_surfStartingParticles.putVar(startingParticles);
//nc_surfZ.putVar(surfZ);
std::cout << "writing energy distribution file " << std::endl;
nc_surfEDist.putVar(&energyDistribution[0]);
//NcVar nc_surfEDistGrid = ncFile1.addVar("gridE",ncDouble,nc_nEnergies);
//nc_surfEDistGrid.putVar(&surfaces->gridE[0]);
//NcVar nc_surfADistGrid = ncFile1.addVar("gridA",ncDouble,nc_nAngles);
//nc_surfADistGrid.putVar(&surfaces->gridA[0]);
ncFile1.close();
#endif
#if PARTICLE_TRACKS > 0

// Write netCDF output for histories
NcFile ncFile_hist("output/history.nc", NcFile::replace);
NcDim nc_nT = ncFile_hist.addDim("nT",nT/subSampleFac);
NcDim nc_nP = ncFile_hist.addDim("nP",nP);
vector<NcDim> dims_hist;
dims_hist.push_back(nc_nP);
dims_hist.push_back(nc_nT);
//NcDim nc_nPnT = ncFile_hist.addDim("nPnT",nP*nT/subSampleFac);
//dims_hist.push_back(nc_nPnT);
NcVar nc_x = ncFile_hist.addVar("x",ncDouble,dims_hist);
NcVar nc_y = ncFile_hist.addVar("y",ncDouble,dims_hist);
NcVar nc_z = ncFile_hist.addVar("z",ncDouble,dims_hist);

NcVar nc_vx = ncFile_hist.addVar("vx",ncDouble,dims_hist);
NcVar nc_vy = ncFile_hist.addVar("vy",ncDouble,dims_hist);
NcVar nc_vz = ncFile_hist.addVar("vz",ncDouble,dims_hist);

NcVar nc_charge = ncFile_hist.addVar("charge",ncDouble,dims_hist);
#if USE_CUDA > 0
nc_x.putVar(&positionHistoryX[0]);
nc_y.putVar(&positionHistoryY[0]);
nc_z.putVar(&positionHistoryZ[0]);

nc_vx.putVar(&velocityHistoryX[0]);
nc_vy.putVar(&velocityHistoryY[0]);
nc_vz.putVar(&velocityHistoryZ[0]);

nc_charge.putVar(&chargeHistory[0]);
#else
nc_x.putVar(positionHistoryX[0]);
nc_y.putVar(positionHistoryY[0]);
nc_z.putVar(positionHistoryZ[0]);

nc_vx.putVar(velocityHistoryX[0]);
nc_vy.putVar(velocityHistoryY[0]);
nc_vz.putVar(velocityHistoryZ[0]);

nc_charge.putVar(chargeHistory[0]);
#endif
ncFile_hist.close();
#endif
#if SPECTROSCOPY > 0
// Write netCDF output for density data
NcFile ncFile("output/spec.nc", NcFile::replace);
NcDim nc_nBins = ncFile.addDim("nBins",nBins+1);
NcDim nc_nR = ncFile.addDim("nR",net_nX);
#if SPECTROSCOPY > 2
NcDim nc_nY = ncFile.addDim("nY",net_nY);
#endif
NcDim nc_nZ = ncFile.addDim("nZ",net_nZ);
vector<NcDim> dims;
dims.push_back(nc_nBins);
dims.push_back(nc_nZ);
#if SPECTROSCOPY > 2
dims.push_back(nc_nY);
#endif
dims.push_back(nc_nR);

NcVar nc_n = ncFile.addVar("n",ncDouble,dims);
NcVar nc_gridR = ncFile.addVar("gridR",ncDouble,nc_nR);
NcVar nc_gridZ = ncFile.addVar("gridZ",ncDouble,nc_nZ);
float *binPointer = &net_Bins[0];
nc_gridR.putVar(&gridX_bins[0]);
nc_gridZ.putVar(&gridZ_bins[0]);
nc_n.putVar(binPointer);
ncFile.close();
#endif
#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif
#if USE_BOOST
    /*
    cpu_times copyToHostTime = timer.elapsed();

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "Copy to host, bin and output time: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';
    std::cout << "Total ODE integration time: " << moveTime*1e-9 << '\n';
    std::cout << "Total geometry checking time: " << geomCheckTime*1e-9 << '\n';
    std::cout << "Total ionization time: " << ionizTime*1e-9 << '\n';
    */
#endif
    for(int i=0;i<100;i++)
{
    std::cout << "particle hitwall and Ez " << particleArray->hitWall[i] << " " << particleArray->test[i] << " "<< particleArray->test0[i] << " " << particleArray->test1[i]<< " "<<
        particleArray->test2[i] << " " << particleArray->test3[i] << std::endl;
}
    for(int i=0;i<100;i++)
{
    std::cout << "particle ionization z and t " << particleArray->firstIonizationZ[i] << " " << particleArray->firstIonizationT[i]  << std::endl;
}
    for(int i=0;i<100;i++)
{
    std::cout << "reflected/sputtered energy " << particleArray->newVelocity[i]   << std::endl;
}
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
    return 0;
    }
