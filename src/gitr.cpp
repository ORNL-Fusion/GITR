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

int main()
{
  //Prepare config files for import
  Config cfg,cfg_geom;
  cfg.readFile("gitrInput.cfg");
  std::cout << "staged input file " << std::endl;
  cfg_geom.readFile("gitrGeometry.cfg");
  std::cout << "staged input files " << std::endl;
  //check binary compatibility with input file
  #if CHECK_COMPATIBILITY>0
    Setting& flags = cfg.lookup("flags");
    
    int check1 = flags["USE_CUDA"];
    std::cout << "check1 "<< typeid(check1).name() << check1 << "use cuda "<< typeid(USE_CUDA).name() << USE_CUDA << std::endl;    
    std::cout << "subtract " << USE_CUDA - check1 << std::endl;
    bool check0 = (USE_CUDA != check1);
    std::cout << "bool " << check0 << std::endl;
    if (USE_CUDA != check1)
    { std::cout << "incompatibility in USE_CUDA between input file and binary" << std::endl;
      exit(0);
    }
  
    int check2 = flags["USEMPI"];      
    if (USEMPI != check2)
    { std::cout << "incompatibility in USEMPI between input file and binary" << std::endl;
      exit(0);
    }

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
    
    int check5 = flags["USERECOMBINATION"];      
    if (USERECOMBINATION != check5)
    { std::cout << "incompatibility in USERECOMBINATION between input file and binary" << std::endl;
      exit(0);
    }
    
    int check6 = flags["USEPERPDIFFUSION"];      
    if (USEPERPDIFFUSION !=check6)
    { std::cout << "incompatibility in USEPERPDIFFUSION between input file and binary" << std::endl;
      exit(0);
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
    int check11 = flags["USEPRESHEATHFIELD"];      
    if (USEPRESHEATHFIELD !=check11)
    {
            std::cout << "incompatibility in USEPRESHEATHFIELD between input file and binary" << std::endl;
            exit(0);
    }
    int check12 = flags["BFIELD_INTERP"];      
    if (BFIELD_INTERP !=check12)
    {
            std::cout << "incompatibility in BFIELD_INTERP between input file and binary" << std::endl;
            exit(0);
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
  float background_Z = cfg.lookup("backgroundPlasmaProfiles.Z");
  float background_amu = cfg.lookup("backgroundPlasmaProfiles.amu");

  //Bfield initialization
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
  float perpDiffusionCoeff = cfg.lookup("backgroundPlasmaProfiles.Diffusion.Dperp");

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
#endif 
float* finalPosX = new float[nP];
float* finalPosY = new float[nP];
float* finalPosZ = new float[nP];
float* finalVx = new float[nP];
float* finalVy = new float[nP];
float* finalVz = new float[nP];
float* transitTime = new float[nP];
#if USE_BOOST
cpu_timer timer;
#endif

std::uniform_real_distribution<float> dist(0,1e6);

#if FIXEDSEEDS == 0
    std::random_device rd;
    std::default_random_engine generator(rd());
#endif
/*
    sim::Array<Particle> particleArray(nParticles);
    for(int i=0;i<nParticles; i++)
    {
        particleArray[i] = particleArray[i];
    }

    for(int i=0;i<nParticles; i++)
    {
        std::cout << particleArray[i].x << " " << particleArray[i].xprevious << std::endl;
        std::cout << particleArray[i].y << " " << particleArray[i].yprevious << std::endl;
        std::cout << particleArray[i].z << " " << particleArray[i].zprevious << std::endl;
        std::cout << particleArray[i].vx << " " << particleArray[i].vy << " " << particleArray[i].vz << std::endl;
        std::cout << particleArray[i].Z << " " << particleArray[i].amu << " " << particleArray[i].charge << std::endl;

    }
    */

      thrust::counting_iterator<std::size_t> particleBegin(0);  
        thrust::counting_iterator<std::size_t> particleEnd(nParticles);
#if PARTICLESEEDS > 0
#if USEIONIZATION > 0
#if FIXEDSEEDS ==1
    float ionization_seeds = cfg.lookup("operators.ionization.seed");
    std::default_random_engine generator(ionization_seeds);
#endif
//#ifdef __CUDACC__
    //sim::Array<Particle> particleArray(nParticles);
    //for(int i=0;i<nParticles; i++)
   // {
   //     std::cout << particleArray[i].x << std::endl;
   // }
    std::cout << "fixed particle seeds " << std::endl;
    sim::Array<float> seeds0(nP);
    std::generate( seeds0.begin(), seeds0.end(), [&]() { return dist(generator); } );
//    thrust::device_vector<float> deviceSeeds0 = seeds0;
    thrust::transform(thrust::device, particleArray->streams.begin(), particleArray->streams.end(),
                    seeds0.begin(), particleArray->streams.begin(), randInit(0) );
/*
#else
    std::vector<float> seeds0(nP);
    std::generate( seeds0.begin(), seeds0.end(), [&]() { return dist(generator); } );
    std::transform(particleArray.begin(), particleArray.end(),
                    seeds0.begin(), particleArray.begin(), randInit(0) );
#endif
*/
#endif

#if USERECOMBINATION > 0
        std::vector<float> seeds1(nP);
        std::generate( seeds1.begin(), seeds1.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds1 = seeds1;
        thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds1.begin(), particleArray.begin(), randInit(1) );
#else
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds1.begin(), particleArray.begin(), randInit(1) );
#endif
#endif

#if USEPERPDIFFUSION > 0
        std::vector<float> seeds2(nP);
        std::generate( seeds2.begin(), seeds2.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds2 = seeds2;
        thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds2.begin(), particleArray.begin(), randInit(2) );
#else
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds2.begin(), particleArray.begin(), randInit(2) );
#endif
#endif

#if USECOULOMBCOLLISIONS > 0
        std::vector<float> seeds3(nP),seeds4(nP),seeds5(nP);
        std::generate( seeds3.begin(), seeds3.end(), [&]() { return dist(generator); } );
    std::generate( seeds4.begin(), seeds4.end(), [&]() { return dist(generator); } );
    std::generate( seeds5.begin(), seeds5.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds3 = seeds3,deviceSeeds4 = seeds4,deviceSeeds5 = seeds5;
        thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds3.begin(), particleArray.begin(), randInit(3) );
    thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds4.begin(), particleArray.begin(), randInit(4) );
        thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds5.begin(), particleArray.begin(), randInit(5) );
#else
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds3.begin(), particleArray.begin(), randInit(3) );
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds4.begin(), particleArray.begin(), randInit(4) );
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds5.begin(), particleArray.begin(), randInit(5) );
#endif
#endif

#if USESURFACEMODEL > 0
        std::vector<float> seeds6(nP);
        std::generate( seeds6.begin(), seeds6.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds6 = seeds6;
        thrust::transform(particleArray.begin(), particleArray.end(),
                    deviceSeeds6.begin(), particleArray.begin(), randInit(6) );
#else
        std::transform(particleArray.begin(), particleArray.end(),
                    seeds6.begin(), particleArray.begin(), randInit(6) );
#endif
#endif

#if __CUDACC__
sim::Array<curandState> state1(7);
#else
sim::Array<std::mt19937> state1(7);
#endif
#else
#if __CUDACC__
sim::Array<curandState> state1(7);
curandInitialize<<<1,1>>>(&state1[0],19);
#else
sim::Array<std::mt19937> state1(7);
std::mt19937 s(348763);
state1[0] = s;
#endif
#endif

    float moveTime = 0.0;
    float geomCheckTime = 0.0;
    float ionizTime = 0.0;

#if USE_BOOST
    cpu_times copyToDeviceTime = timer.elapsed();
    std::cout << "Initialize rand state and copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';
#endif
    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::duration<float> fsec;
    auto start_clock = Time::now();
    std::cout << "Starting main loop" << std::endl;
//Main time loop
    for(int tt=0; tt< nT; tt++)
    {
//#ifdef __CUDACC__
        thrust::for_each(thrust::device, particleBegin,particleEnd, 
                move_boris(particleArray,dt,boundaries.data(), nLines,
                    nR_Bfield,nZ_Bfield, bfieldGridr.data(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front(),
                    nR_PreSheathEfield,nZ_PreSheathEfield,
                    &preSheathEGridr.front(),&preSheathEGridz.front(),
                    &PSEr.front(),&PSEz.front(),&PSEt.front(),
                        nR_closeGeom_sheath,nZ_closeGeom_sheath,n_closeGeomElements_sheath,
                        &closeGeomGridr_sheath.front(),&closeGeomGridz_sheath.front(),
                        &closeGeom_sheath.front()) );
        
        
        //try {
            thrust::for_each(thrust::device, particleBegin,particleEnd,
                    geometry_check(particleArray,boundaryModArray,nLines,&boundaries[0],dt,tt,
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
                    &gridDensity_Ionization.front(), &rateCoeff_Ionization.front(),tt));
#endif
#if USERECOMBINATION > 0
        thrust::for_each(particleArray.begin(), particleArray.end(),
                recombine(dt) );
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
        thrust::for_each(particleArray.begin(), particleArray.end(), 
                reflection(dt,nLines,BoundaryDevicePointer) );
#endif        

#if PARTICLE_TRACKS >0
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
    std::for_each(particleArray.begin(), particleArray.end(),
            ionize(dt,&state1.front(),
                nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),
                nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                nTemperaturesIonize, 
                nDensitiesIonize, &gridTemperature_Ionization.front(),
                &gridDensity_Ionization.front(), &rateCoeff_Ionization.front() ) );
#if USE_BOOST
cpu_times ionizTime1 = timer.elapsed();
ionizTime = ionizTime + (ionizTime1.wall - ionizTime0.wall);
#endif
#endif
#if USERECOMBINATION > 0
    std::for_each(particleArray.begin(), particleArray.end(), 
            recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
        std::for_each(particleArray.begin(), particleArray.end(), 
                crossFieldDiffusion(dt,perpDiffusionCoeff,
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front()));
        
        std::for_each(particleArray.begin(), particleArray.end(), 
                geometry_check(nLines,boundaries.data(),dt,tt) );
#endif
#if USECOULOMBCOLLISIONS > 0
        std::for_each(particleArray.begin(), particleArray.end(), 
                coulombCollisions(dt,
                    nR_flowV,nZ_flowV,&flowVGridr.front(),&flowVGridz.front(),
                    &flowVr.front(),&flowVz.front(),&flowVt.front(),
                    nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),
                    nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                    background_Z,background_amu,
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),
                    &bfieldGridz.front(),&br.front(),&bz.front(),&bt.front()));
#endif
#if USETHERMALFORCE > 0
        std::for_each(particleArray.begin(), particleArray.end(), 
                thermalForce(dt,background_amu,
                    nR_gradT,nZ_gradT,&gradTGridr.front(),&gradTGridz.front(),
                    &gradTiR.front(),&gradTiZ.front(),&gradTiT.front(),
                    &gradTeR.front(),&gradTeZ.front(),&gradTeT.front(),
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front()));
#endif
#if USESURFACEMODEL > 0
        std::for_each(particleArray.begin(), particleArray.end(), 
                reflection(dt,nLines,boundaries.data()) );
#endif        
#if PARTICLE_TRACKS >0
if (tt % subSampleFac == 0)  
{    
        for(int i=0;i<nP;i++)
        {
            positionHistoryX[i][tt/subSampleFac] = particleArray[i].xprevious;
            positionHistoryY[i][tt/subSampleFac] = particleArray[i].yprevious;
            positionHistoryZ[i][tt/subSampleFac] = particleArray[i].zprevious;
            velocityHistoryX[i][tt/subSampleFac] = particleArray[i].vx;
            velocityHistoryY[i][tt/subSampleFac] = particleArray[i].vy;
            velocityHistoryZ[i][tt/subSampleFac] = particleArray[i].vz;
        }
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
#if USE3DTETGEOM > 0
    float meanTransitTime0 = 0.0;
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

NcVar nc_surfImpacts = ncFile1.addVar("impacts",ncDouble,dims1);
nc_surfImpacts.putVar(impacts);
#if PARTICLE_TRACKS > 0

// Write netCDF output for histories
NcFile ncFile_hist("history.nc", NcFile::replace);
NcDim nc_nT = ncFile_hist.addDim("nT",nT/subSampleFac);
NcDim nc_nP = ncFile_hist.addDim("nP",nP);
vector<NcDim> dims_hist;
dims_hist.push_back(nc_nP);
dims_hist.push_back(nc_nT);

NcVar nc_x = ncFile_hist.addVar("x",ncDouble,dims_hist);
NcVar nc_y = ncFile_hist.addVar("y",ncDouble,dims_hist);
NcVar nc_z = ncFile_hist.addVar("z",ncDouble,dims_hist);

NcVar nc_vx = ncFile_hist.addVar("vx",ncDouble,dims_hist);
NcVar nc_vy = ncFile_hist.addVar("vy",ncDouble,dims_hist);
NcVar nc_vz = ncFile_hist.addVar("vz",ncDouble,dims_hist);

NcVar nc_charge = ncFile_hist.addVar("charge",ncDouble,dims_hist);

nc_x.putVar(positionHistoryX[0]);
nc_y.putVar(positionHistoryY[0]);
nc_z.putVar(positionHistoryZ[0]);

nc_vx.putVar(velocityHistoryX[0]);
nc_vy.putVar(velocityHistoryY[0]);
nc_vz.putVar(velocityHistoryZ[0]);

nc_charge.putVar(chargeHistory[0]);
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

#ifdef __CUDACC__
    cudaError_t err = cudaDeviceReset();
//cudaProfilerStop();
#endif
    return 0;
    }
