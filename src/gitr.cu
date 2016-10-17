#include <iostream>
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
#include "Particle.h"
#include "Boundary.h"
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

#ifdef __CUDACC__
    #include <thrust/copy.h>
    #include <thrust/host_vector.h>
    #include <thrust/device_vector.h>
    #include <curand.h>
    #include <curand_kernel.h>
    #include <thrust/sequence.h>
    #include <thrust/transform.h>
    #include <thrust/functional.h>
    #include <thrust/device_ptr.h>
#endif

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
cfg_geom.readFile("gitrGeometry.cfg");


// Background species info
float background_Z = cfg.lookup("backgroundPlasmaProfiles.Z");
float background_amu = cfg.lookup("backgroundPlasmaProfiles.amu");

//Bfield initialization
#if BFIELD_INTERP == 0
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    std::vector<float> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    std::vector<float> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);
    br[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.br");
    bz[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bz");
    bt[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bt");
#elif BFIELD_INTERP == 2
    int nR_Bfield;
    int nZ_Bfield;
    
    int b1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNzString"),nR_Bfield,nZ_Bfield);
    
    std::vector<float> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    std::vector<float> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);
    
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

#ifdef __CUDACC__
    thrust::device_vector<float> deviceBfieldGridRVector = bfieldGridr;
    thrust::device_vector<float> deviceBfieldGridZVector = bfieldGridz;
    thrust::device_vector<float> deviceBfieldRVector = br;
    thrust::device_vector<float> deviceBfieldZVector = bz;
    thrust::device_vector<float> deviceBfieldTVector = bt;
    
    float * BfieldGridRDevicePointer = thrust::raw_pointer_cast(deviceBfieldGridRVector.data());
    float * BfieldGridZDevicePointer = thrust::raw_pointer_cast(deviceBfieldGridZVector.data());
    float * BfieldRDevicePointer = thrust::raw_pointer_cast(deviceBfieldRVector.data());
    float * BfieldZDevicePointer = thrust::raw_pointer_cast(deviceBfieldZVector.data());
    float * BfieldTDevicePointer = thrust::raw_pointer_cast(deviceBfieldTVector.data());
#endif

//Background Plasma Temperature Initialization    
#if TEMP_INTERP == 0
    int nR_Temp = 1;
    int nZ_Temp = 1;
    std::vector<float> TempGridr(nR_Temp), TempGridz(nZ_Temp);
    std::vector<float> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);
    ti[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.ti");
    te[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.te");
#elif TEMP_INTERP == 2
    int nR_Temp;
    int nZ_Temp;
    
    int t1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNzString"),nR_Temp,nZ_Temp);
    std::vector<float> TempGridr(nR_Temp), TempGridz(nZ_Temp);
    std::vector<float> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);
    
    int t2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridRString"), TempGridr);
    
    int t3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridZString"), TempGridz);
    std::cout << "temperature import" << nZ_Temp << nR_Temp << std::endl;
    int t4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.IonTempString"), ti);
    
    int t5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.ElectronTempString"), te);
#endif
std::string outnameTi = "ti.m";
std::string outnameTe = "te.m";
OUTPUT2d(profiles_folder,outnameTi, nR_Temp, nZ_Temp, &ti.front());
OUTPUT2d(profiles_folder,outnameTe, nR_Temp, nZ_Temp, &te.front());

#ifdef __CUDACC__
    thrust::device_vector<float> deviceTempGridRVector = TempGridr;
    thrust::device_vector<float> deviceTempGridZVector = TempGridz;
    thrust::device_vector<float> deviceTiVector = ti;
    thrust::device_vector<float> deviceTeVector = te;
    float * TempGridRDevicePointer = thrust::raw_pointer_cast(deviceTempGridRVector.data());
    float * TempGridZDevicePointer = thrust::raw_pointer_cast(deviceTempGridZVector.data());
    float * TiDevicePointer = thrust::raw_pointer_cast(deviceTiVector.data());
    float * TeDevicePointer = thrust::raw_pointer_cast(deviceTeVector.data());
#endif

//Background Plasma Density Initialization
#if DENSITY_INTERP == 0
    int nR_Dens = 1;
    int nZ_Dens = 1;
    std::vector<float> DensGridr(nR_Dens), DensGridz(nZ_Dens);
    std::vector<float> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
    ni[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ni");
    ne[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ne");
#elif DENSITY_INTERP == 2
    int nR_Dens;
    int nZ_Dens;
    
    int n1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Density.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Density.gridNzString"),nR_Dens,nZ_Dens);
    
    std::vector<float> DensGridr(nR_Dens), DensGridz(nZ_Dens);
    std::vector<float> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
    
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

#ifdef __CUDACC__
    thrust::device_vector<float> deviceDensGridRVector = DensGridr;
    thrust::device_vector<float> deviceDensGridZVector = DensGridz;
    thrust::device_vector<float> deviceNiVector = ni;
    thrust::device_vector<float> deviceNeVector = ne;
    float * DensGridRDevicePointer = thrust::raw_pointer_cast(deviceDensGridRVector.data());
    float * DensGridZDevicePointer = thrust::raw_pointer_cast(deviceDensGridZVector.data());
    float * NiDevicePointer = thrust::raw_pointer_cast(deviceNiVector.data());
    float * NeDevicePointer = thrust::raw_pointer_cast(deviceNeVector.data());
#endif

//Connection length initialization
#if (FLOWV_INTERP == 1 || FLOWV_INTERP == 4)
    int nR_Lc;
    int nZ_Lc;
    
    int l1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.gridNzString"),nR_Lc,nZ_Lc);
    std::vector<float> LcGridr(nR_Lc), LcGridz(nZ_Lc);
    std::vector<float> Lc(nR_Lc*nZ_Lc), s(nR_Lc*nZ_Lc);
    
    int l2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.gridRString"), LcGridr);
    
    int l3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.gridZString"), LcGridz);
    
    int l4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.LcString"), Lc);
    
    //int l5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.fileString"),
            //cfg.lookup("backgroundPlasmaProfiles.ConnectionLength.SString"), s);
#endif
std::string outnameLc = "Lc.m";
std::string outnameS = "s.m";
OUTPUT2d(profiles_folder,outnameLc, nR_Lc, nZ_Lc, &Lc.front());
OUTPUT2d(profiles_folder,outnameS, nR_Lc, nZ_Lc, &s.front());

#ifdef __CUDACC__
    thrust::device_vector<float> deviceLcGridRVector = LcGridr;
    thrust::device_vector<float> deviceLcGridZVector = LcGridz;
    thrust::device_vector<float> deviceLcVector = Lc;
    thrust::device_vector<float> deviceSVector = s;
    float * LcGridRDevicePointer = thrust::raw_pointer_cast(deviceLcGridRVector.data());
    float * LcGridZDevicePointer = thrust::raw_pointer_cast(deviceLcGridZVector.data());
    float * LcDevicePointer = thrust::raw_pointer_cast(deviceLcVector.data());
    float * SDevicePointer = thrust::raw_pointer_cast(deviceSVector.data());
#endif

//Background Plasma flow velocity initialization    
#if FLOWV_INTERP == 0
    int nR_flowV = 1;
    int nZ_flowV = 1;
    std::vector<float> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
    std::vector<float> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),flowVt(nR_flowV*nZ_flowV);
    flowVr[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVr");
    flowVz[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVz");
#elif (FLOWV_INTERP == 2 || FLOWV_INTERP == 4)
    int nR_flowV;
    int nZ_flowV;
    
    int f1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNzString"),nR_flowV,nZ_flowV);
    
    std::vector<float> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
    std::vector<float> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),flowVt(nR_flowV*nZ_flowV);
    
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

#ifdef __CUDACC__
    thrust::device_vector<float> deviceFlowVGridRVector = flowVGridr;
    thrust::device_vector<float> deviceFlowVGridZVector = flowVGridz;
    thrust::device_vector<float> deviceFlowVrVector = flowVr;
    thrust::device_vector<float> deviceFlowVzVector = flowVz;
    thrust::device_vector<float> deviceFlowVtVector = flowVt;
    float * FlowVGridRDevicePointer = thrust::raw_pointer_cast(deviceFlowVGridRVector.data());
    float * FlowVGridZDevicePointer = thrust::raw_pointer_cast(deviceFlowVGridZVector.data());
    float * FlowVrDevicePointer = thrust::raw_pointer_cast(deviceFlowVrVector.data());
    float * FlowVzDevicePointer = thrust::raw_pointer_cast(deviceFlowVzVector.data());
    float * FlowVtDevicePointer = thrust::raw_pointer_cast(deviceFlowVtVector.data());
#endif

//Background plasma temperature gradient field intitialization    
#if GRADT_INTERP == 0
    int nR_gradT = 1;
    int nZ_gradT = 1;
    std::vector<float> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
    std::vector<float> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),gradTeT(nR_gradT*nZ_gradT,0.0),
        gradTiR(nR_gradT*nZ_gradT), gradTiZ(nR_gradT*nZ_gradT),gradTiT(nR_gradT*nZ_gradT,0.0);
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
    
    std::vector<float> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
    std::vector<float> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),gradTeT(nR_gradT*nZ_gradT,0.0),
        gradTiR(nR_gradT*nZ_gradT), gradTiZ(nR_gradT*nZ_gradT),gradTiT(nR_gradT*nZ_gradT,0.0);
    
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

#ifdef __CUDACC__
    thrust::device_vector<float> deviceGradTGridRVector = gradTGridr;
    thrust::device_vector<float> deviceGradTGridZVector = gradTGridz;
    thrust::device_vector<float> deviceGradTiRVector = gradTiR;
    thrust::device_vector<float> deviceGradTiZVector = gradTiZ;
    thrust::device_vector<float> deviceGradTiTVector = gradTiT;
    thrust::device_vector<float> deviceGradTeRVector = gradTeR;
    thrust::device_vector<float> deviceGradTeZVector = gradTeZ;
    thrust::device_vector<float> deviceGradTeTVector = gradTeT;
    float * GradTGridRDevicePointer = thrust::raw_pointer_cast(deviceGradTGridRVector.data());
    float * GradTGridZDevicePointer = thrust::raw_pointer_cast(deviceGradTGridZVector.data());
    float * GradTiRDevicePointer = thrust::raw_pointer_cast(deviceGradTiRVector.data());
    float * GradTiZDevicePointer = thrust::raw_pointer_cast(deviceGradTiZVector.data());
    float * GradTiTDevicePointer = thrust::raw_pointer_cast(deviceGradTiTVector.data());
    float * GradTeRDevicePointer = thrust::raw_pointer_cast(deviceGradTeRVector.data());
    float * GradTeZDevicePointer = thrust::raw_pointer_cast(deviceGradTeZVector.data());
    float * GradTeTDevicePointer = thrust::raw_pointer_cast(deviceGradTeTVector.data());
#endif

//Initialization of ionization and recombination coefficients    
int nCS_Ionize, nCS_Recombine;
int i0 = read_profileNs(cfg.lookup("impurityParticleSource.ionization.fileString"),
            cfg.lookup("impurityParticleSource.ionization.nChargeStateString"),
            cfg.lookup("impurityParticleSource.recombination.nChargeStateString"),
            nCS_Ionize, nCS_Recombine);
int nTemperaturesIonize;
int nDensitiesIonize;
int i1 = read_profileNs(cfg.lookup("impurityParticleSource.ionization.fileString"),
            cfg.lookup("impurityParticleSource.ionization.DensGridString"),
            cfg.lookup("impurityParticleSource.ionization.TempGridString"),
            nDensitiesIonize,nTemperaturesIonize);

std::vector<float> rateCoeff_Ionization(nCS_Ionize*nTemperaturesIonize*nDensitiesIonize);
std::vector<float> gridTemperature_Ionization(nTemperaturesIonize),
                        gridDensity_Ionization(nDensitiesIonize);

int i2 = read_profiles(cfg.lookup("impurityParticleSource.ionization.fileString"),
        nTemperaturesIonize,nDensitiesIonize,
        cfg.lookup("impurityParticleSource.ionization.TempGridVarName"), 
        gridTemperature_Ionization,cfg.lookup("impurityParticleSource.ionization.DensGridVarName"),
        gridDensity_Ionization,
        cfg.lookup("impurityParticleSource.ionization.CoeffVarName"),
        rateCoeff_Ionization);
   
int nTemperaturesRecombine;
int nDensitiesRecombine;
int i3 = read_profileNs(cfg.lookup("impurityParticleSource.recombination.fileString"),
            cfg.lookup("impurityParticleSource.recombination.DensGridString"),
            cfg.lookup("impurityParticleSource.recombination.TempGridString"),
            nDensitiesRecombine,nTemperaturesRecombine);

std::vector<float> rateCoeff_Recombination(nCS_Recombine*nTemperaturesRecombine*nDensitiesRecombine);
std::vector<float> gridTemperature_Recombination(nTemperaturesRecombine),
                    gridDensity_Recombination(nDensitiesRecombine);

int i4 = read_profiles(cfg.lookup("impurityParticleSource.recombination.fileString"),
             nTemperaturesRecombine,nDensitiesRecombine,
             cfg.lookup("impurityParticleSource.recombination.TempGridVarName"), 
             gridTemperature_Recombination,cfg.lookup("impurityParticleSource.recombination.DensGridVarName"),
             gridDensity_Recombination,
             cfg.lookup("impurityParticleSource.recombination.CoeffVarName"),
             rateCoeff_Recombination);

#ifdef __CUDACC__
    thrust::device_vector<float> deviceTemperatureGridIonizationVector = gridTemperature_Ionization;
    thrust::device_vector<float> deviceDensityGridIonizationVector = gridDensity_Ionization;
    thrust::device_vector<float> deviceTemperatureGridRecombinationVector = gridTemperature_Recombination;
    thrust::device_vector<float> deviceDensityGridRecombinationVector = gridDensity_Recombination;
    thrust::device_vector<float> deviceRateCoeffIonizationVector = rateCoeff_Ionization;
    thrust::device_vector<float> deviceRateCoeffRecombinationVector = rateCoeff_Recombination;
    float * IonizationTemperatureGridDevicePointer = thrust::raw_pointer_cast(deviceTemperatureGridIonizationVector.data());
    float * RecombinationTemperatureGridDevicePointer = thrust::raw_pointer_cast(deviceTemperatureGridRecombinationVector.data());
    float * IonizationDensityGridDevicePointer = thrust::raw_pointer_cast(deviceDensityGridIonizationVector.data());
    float * RecombinationDensityGridDevicePointer = thrust::raw_pointer_cast(deviceDensityGridRecombinationVector.data());
    float * IonizationRateCoeffDevicePointer = thrust::raw_pointer_cast(deviceRateCoeffIonizationVector.data());
    float * RecombinationRateCoeffDevicePointer = thrust::raw_pointer_cast(deviceRateCoeffRecombinationVector.data());
#endif

//Geometry Definition
Setting& geom = cfg_geom.lookup("geom");
int nLines = geom["x1"].getLength();
int nMaterials = geom["nMaterials"];
std::cout << "Number of Geometric Objects Loaded: " << nLines << std::endl;

std::vector<Boundary> hostBoundaryVector(nLines+1);

std::string geom_outname = "geom.m";
std::string geom_folder = "geometry";
ofstream outfile;

#if USE_BOOST
//Output
boost::filesystem::path dir(geom_folder);

if(!(boost::filesystem::exists(dir)))
{
 // std::cout<<"Doesn't Exists"<<std::endl;
      if (boost::filesystem::create_directory(dir))
      {
      //std::cout << " Successfully Created " << std::endl;
      }
}
#endif
std::string full_path = geom_folder + "/" + geom_outname;
outfile.open (full_path );
for(int i=0; i<nMaterials; i++)
{
}
for(int i=0 ; i<nLines ; i++)
{
     hostBoundaryVector[i].x1 = geom["x1"][i];
     hostBoundaryVector[i].z1 = geom["z1"][i];
     hostBoundaryVector[i].x2 = geom["x2"][i];
     hostBoundaryVector[i].z2 = geom["z2"][i];
     hostBoundaryVector[i].Z = geom["Z"][i];
     hostBoundaryVector[i].slope_dzdx = geom["slope"][i];
     hostBoundaryVector[i].intercept_z = geom["intercept"][i];
     hostBoundaryVector[i].length = geom["length"][i];

     outfile << "geom(" << i+1 << ",:) = ["<<hostBoundaryVector[i].x1 << ", " <<hostBoundaryVector[i].z1 << ", " <<
        hostBoundaryVector[i].x2 << ", " << hostBoundaryVector[i].z2 << ", " <<
        hostBoundaryVector[i].slope_dzdx << ", " << hostBoundaryVector[i].intercept_z << ", " <<
        hostBoundaryVector[i].length << ", " << hostBoundaryVector[i].Z << "];" << std::endl;
}   

outfile.close();

hostBoundaryVector[nLines].Z = geom["Z"][nLines];
hostBoundaryVector[nLines].y1 = geom["y1"];
hostBoundaryVector[nLines].y2 = geom["y2"];
hostBoundaryVector[nLines].periodic = geom["periodic"];
std::cout << "Starting Boundary Init..." << std::endl;

//Applying background values at material boundaries
std::for_each(hostBoundaryVector.begin(), hostBoundaryVector.end()-1,
            boundary_init(background_Z,background_amu,
            nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ni.front(),
            nR_Bfield,nZ_Bfield,&bfieldGridr.front(),
            &bfieldGridz.front(),&br.front(),&bz.front(), &bt.front(),
            nR_Temp,nZ_Temp,&TempGridr.front(),
            &TempGridz.front(),&ti.front() ));

#ifdef __CUDACC__
    thrust::device_vector<Boundary> deviceBoundaryVector = hostBoundaryVector;
    Boundary * BoundaryDevicePointer = thrust::raw_pointer_cast(deviceBoundaryVector.data());
#endif
    std::cout << "Completed Boundary Init " << std::endl;
//Efield
#if USEPRESHEATHEFIELD > 0    
#if PRESHEATH_INTERP == 0
    int nR_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;
    std::vector<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
    std::vector<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
        PSEt(nR_PreSheathEfield*nZ_PreSheathEfield);
    PSEr[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Er");
    PSEz[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Ez");
    PSEt[0] = cfg.lookup("backgroundPlasmaProfiles.Efield.Et");
#elif PRESHEATH_INTERP == 2
    int nR_PreSheathEfield;
    int nZ_PreSheathEfield;
    
    int e1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Efield.fileString"),
                cfg.lookup("backgroundPlasmaProfiles.Efield.gridNrString"),
                cfg.lookup("backgroundPlasmaProfiles.Efield.gridNzString"),nR_PreSheathEfield,nZ_PreSheathEfield);
    
    std::vector<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
    std::vector<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
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
    std::vector<float> preSheathEGridr(nR_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);
    std::vector<float> PSEr(nR_PreSheathEfield*nZ_PreSheathEfield), PSEz(nR_PreSheathEfield*nZ_PreSheathEfield),
        PSEt(nR_PreSheathEfield*nZ_PreSheathEfield);
#endif
#ifdef __CUDACC__
    thrust::device_vector<float> deviceEfieldGridRVector = preSheathEGridr;
    thrust::device_vector<float> deviceEfieldGridZVector = preSheathEGridz;
    thrust::device_vector<float> deviceEfieldRVector = PSEr;
    thrust::device_vector<float> deviceEfieldZVector = PSEz;
    thrust::device_vector<float> deviceEfieldTVector = PSEt;
    
    float * EfieldGridRDevicePointer = thrust::raw_pointer_cast(deviceEfieldGridRVector.data());
    float * EfieldGridZDevicePointer = thrust::raw_pointer_cast(deviceEfieldGridZVector.data());
    float * EfieldRDevicePointer = thrust::raw_pointer_cast(deviceEfieldRVector.data());
    float * EfieldZDevicePointer = thrust::raw_pointer_cast(deviceEfieldZVector.data());
    float * EfieldTDevicePointer = thrust::raw_pointer_cast(deviceEfieldTVector.data());
#endif
    std::cout << "Completed presheath Efield Init " << std::endl;
std::vector<float> Efieldr(nR_Bfield*nZ_Bfield), Efieldz(nR_Bfield*nZ_Bfield),Efieldt(nR_Bfield*nZ_Bfield),
                        minDist(nR_Bfield*nZ_Bfield);

#if USESHEATHEFIELD > 0
#if EFIELD_INTERP == 1
float thisE[3] = {0.0,0.0,0.0};
for(int i=0;i<nR_Bfield;i++)
{
    for(int j=0;j<nZ_Bfield;j++)
    {
        minDist[(nR_Bfield - 1 -i)*nZ_Bfield+(nZ_Bfield -1-j)] = getE ( bfieldGridr[i], 0.0, bfieldGridz[j],
                                                                  thisE, hostBoundaryVector.data(),nLines );
        Efieldr[i*nZ_Bfield+j] = thisE[0];
        Efieldz[i*nZ_Bfield+j] = thisE[2];
        Efieldt[i*nZ_Bfield+j] = thisE[1];
    }
}
#endif
#endif
std::string outnameEfieldR = "EfieldR.m";
std::string outnameEfieldZ = "EfieldZ.m";
std::string outnameEfieldT = "EfieldT.m";
std::string outnameMinDist = "DistToSurface.m";
OUTPUT2d(profiles_folder,outnameEfieldR, nR_Bfield, nZ_Bfield, &Efieldr.front());
OUTPUT2d(profiles_folder,outnameEfieldZ, nR_Bfield, nZ_Bfield, &Efieldz.front());
OUTPUT2d(profiles_folder,outnameEfieldT, nR_Bfield, nZ_Bfield, &Efieldt.front());
OUTPUT2d(profiles_folder,outnameMinDist, nR_Bfield, nZ_Bfield, &minDist.front());
    std::cout << "Completed Debye Sheath Init " << std::endl;

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

float x = cfg.lookup("impurityParticleSource.initialConditions.x_start");
float y = cfg.lookup("impurityParticleSource.initialConditions.y_start");
float z = cfg.lookup("impurityParticleSource.initialConditions.z_start");

float Ex = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
float Ey = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
float Ez = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");

float amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
float Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");
float charge = cfg.lookup("impurityParticleSource.initialConditions.charge");
/*
float **SurfaceBins;
float **SurfaceBinsCharge;
float **SurfaceBinsEnergy;
float **SurfaceBinsErosion;

SurfaceBins = new float*[nY];
SurfaceBinsCharge = new float*[nY];
SurfaceBinsEnergy = new float*[nY];
SurfaceBinsErosion = new float*[nY];

SurfaceBins[0] = new float[nY*nZ];
SurfaceBinsCharge[0] = new float[nY*nZ];
SurfaceBinsEnergy[0] = new float[nY*nZ];
SurfaceBinsErosion[0] = new float[nY*nZ];
        
for(int i=0 ; i<nY ; i++)
{
    SurfaceBins[i] = &SurfaceBins[0][i*nZ];
    SurfaceBinsCharge[i] = &SurfaceBinsCharge[0][i*nZ];
    SurfaceBinsEnergy[i] = &SurfaceBinsEnergy[0][i*nZ];
    SurfaceBinsErosion[i] = &SurfaceBinsErosion[0][i*nZ];               
    for(int j=0 ; j<nZ ; j++)
    {
        SurfaceBins[i][j] = 0;
        SurfaceBinsCharge[i][j] = 0;
        SurfaceBinsEnergy[i][j] = 0;
        SurfaceBinsErosion[i][j] = 0;
    }
}
*/    
float dt;
float nPtsPerGyroOrbit = cfg.lookup("timeStep.nPtsPerGyroOrbit");
dt = 2.4e-7/100.0;

const int nP = cfg.lookup("impurityParticleSource.nP");
cout << "Number of particles: " << nP << endl;              
long nParticles = nP;
int nT = cfg.lookup("timeStep.nT");
cout << "Number of time steps: " << nT << " With dt = " << dt << endl; 

//    int surfaceIndexY;
//    int surfaceIndexZ;
#if PARTICLE_SOURCE == 0
    Particle p1(x,y,z,Ex,Ey,Ez,Z,amu,charge);
#ifdef __CUDACC__
    thrust::host_vector<Particle> hostCudaParticleVector(nParticles,p1);
#else
    std::vector<Particle> hostCudaParticleVector(nParticles,p1);
#endif
#elif PARTICLE_SOURCE == 1
    float impurity_Z = cfg.lookup("impurityParticleSource.Z");
    int nImpurityBoundaries = 0;
    for (int i=0; i<nLines;i++)
    {
        if(hostBoundaryVector[i].Z == impurity_Z)
        {
            nImpurityBoundaries++;
        }
    }
    std::cout << "n Impurity Boundaries to launch from " << nImpurityBoundaries << std::endl;
    std::vector<int> boundaryIndex_ImpurityLaunch(nImpurityBoundaries);

    int count = 0;
    for (int i=0; i<nLines;i++)
    {
        if(hostBoundaryVector[i].Z == impurity_Z)
        {
            boundaryIndex_ImpurityLaunch[count] = i;
            count++;
            std::cout << "Boundary indices " << i << std::endl;
        }
    }
    
    int impuritiesPerBoundary = nP/nImpurityBoundaries;
#ifdef __CUDACC__
      thrust::host_vector<Particle> hostCudaParticleVector(nParticles);
#else
      std::vector<Particle> hostCudaParticleVector(nParticles);
#endif
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
            x = hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].x1 + 
                hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].length*rand0;//1.4290;
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
            x = hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].x1 + hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].length*rand0;
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
            
        Particle p1(x,0.0,z,Ex,Ey,Ez,74,184.0,charge);
            hostCudaParticleVector[i*impuritiesPerBoundary + j] = p1;
            //std::cout << " E0 " << E0 << std::endl;
            //std::cout << "vy " << hostCudaParticleVector[i*impuritiesPerBoundary + j].vy << " " << Ey << std::endl;
            //std::cout << "vx " << hostCudaParticleVector[i*impuritiesPerBoundary + j].vx << " " << Ex << std::endl;
            //std::cout << "vz " << hostCudaParticleVector[i*impuritiesPerBoundary + j].vz << " " << Ez << std::endl;
        }
    }
#endif


#if GEOM_TRACE > 0       
    std::uniform_real_distribution<float> dist2(0,1);
    std::random_device rd2;
    std::cout << "Randomizing velocities to trace geometry. " << std::endl;

    for (int i=0 ; i<nParticles ; i++)
    {   float theta = dist2(rd2)*2*3.1415;
        float phi = dist2(rd2)*3.1415;
        float mag = 2e3;
        hostCudaParticleVector[i].vx = mag*cos(theta)*sin(phi);
        hostCudaParticleVector[i].vy = mag*sin(theta)*sin(phi);
        hostCudaParticleVector[i].vz = mag*cos(phi);
    }
#endif

#if PARTICLE_TRACKS > 0
    int subSampleFac = 10;
    float **positionHistoryX;
    float **positionHistoryY;
    float **positionHistoryZ;
    float **velocityHistoryX;
    float **velocityHistoryY;
    float **velocityHistoryZ;
    positionHistoryX = new float* [nP];
    positionHistoryY = new float* [nP];
    positionHistoryZ = new float* [nP];
    velocityHistoryX = new float* [nP];
    velocityHistoryY = new float* [nP];
    velocityHistoryZ = new float* [nP];
    positionHistoryX[0] = new float [nT*nP/subSampleFac];
    positionHistoryY[0] = new float [nT*nP/subSampleFac];
    positionHistoryZ[0] = new float [nT*nP/subSampleFac];
    velocityHistoryX[0] = new float [nT*nP/subSampleFac];
    velocityHistoryY[0] = new float [nT*nP/subSampleFac];
    velocityHistoryZ[0] = new float [nT*nP/subSampleFac];
    for(int i=0 ; i<nP ; i++)
    {
        positionHistoryX[i] = &positionHistoryX[0][i*nT/subSampleFac];
        positionHistoryY[i] = &positionHistoryY[0][i*nT/subSampleFac];
        positionHistoryZ[i] = &positionHistoryZ[0][i*nT/subSampleFac];
        velocityHistoryX[i] = &velocityHistoryX[0][i*nT/subSampleFac];
        velocityHistoryY[i] = &velocityHistoryY[0][i*nT/subSampleFac];
        velocityHistoryZ[i] = &velocityHistoryZ[0][i*nT/subSampleFac];
        for(int j=0 ; j<nT/subSampleFac ; j++)
        {
            positionHistoryX[i][j] = 0;
            positionHistoryY[i][j] = 0;
            positionHistoryZ[i][j] = 0;
            velocityHistoryX[i][j] = 0;
            velocityHistoryY[i][j] = 0;
            velocityHistoryZ[i][j] = 0;
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
#ifdef __CUDACC__
    thrust::device_vector<Particle> deviceCudaParticleVector = hostCudaParticleVector;
#endif

std::uniform_real_distribution<float> dist(0,1e6);

#if FIXEDSEEDS == 0
    std::random_device rd;
    std::default_random_engine generator(rd());
#endif

#if USEIONIZATION > 0
#if FIXEDSEEDS ==1
    float ionization_seeds = cfg.lookup("operators.ionization.seed");
    std::default_random_engine generator(ionization_seeds);
#endif
    std::vector<float> seeds0(nP);
    std::generate( seeds0.begin(), seeds0.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
    thrust::device_vector<float> deviceSeeds0 = seeds0;
    thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds0.begin(), deviceCudaParticleVector.begin(), randInit(0) );
#else
    std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds0.begin(), hostCudaParticleVector.begin(), randInit(0) );
#endif
#endif

#if USERECOMBINATION > 0
        std::vector<float> seeds1(nP);
        std::generate( seeds1.begin(), seeds1.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds1 = seeds1;
        thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds1.begin(), deviceCudaParticleVector.begin(), randInit(1) );
#else
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds1.begin(), hostCudaParticleVector.begin(), randInit(1) );
#endif
#endif

#if USEPERPDIFFUSION > 0
        std::vector<float> seeds2(nP);
        std::generate( seeds2.begin(), seeds2.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds2 = seeds2;
        thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds2.begin(), deviceCudaParticleVector.begin(), randInit(2) );
#else
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds2.begin(), hostCudaParticleVector.begin(), randInit(2) );
#endif
#endif

#if USECOULOMBCOLLISIONS > 0
        std::vector<float> seeds3(nP),seeds4(nP),seeds5(nP);
        std::generate( seeds3.begin(), seeds3.end(), [&]() { return dist(generator); } );
    std::generate( seeds4.begin(), seeds4.end(), [&]() { return dist(generator); } );
    std::generate( seeds5.begin(), seeds5.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds3 = seeds3,deviceSeeds4 = seeds4,deviceSeeds5 = seeds5;
        thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds3.begin(), deviceCudaParticleVector.begin(), randInit(3) );
    thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds4.begin(), deviceCudaParticleVector.begin(), randInit(4) );
        thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds5.begin(), deviceCudaParticleVector.begin(), randInit(5) );
#else
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds3.begin(), hostCudaParticleVector.begin(), randInit(3) );
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds4.begin(), hostCudaParticleVector.begin(), randInit(4) );
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds5.begin(), hostCudaParticleVector.begin(), randInit(5) );
#endif
#endif

#if USESURFACEMODEL > 0
        std::vector<float> seeds6(nP);
        std::generate( seeds6.begin(), seeds6.end(), [&]() { return dist(generator); } );
#ifdef __CUDACC__
        thrust::device_vector<float> deviceSeeds6 = seeds6;
        thrust::transform(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    deviceSeeds6.begin(), deviceCudaParticleVector.begin(), randInit(6) );
#else
        std::transform(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                    seeds6.begin(), hostCudaParticleVector.begin(), randInit(6) );
#endif
#endif
    float moveTime = 0.0;
    float geomCheckTime = 0.0;
    float ionizTime = 0.0;
#if USE_BOOST
    cpu_times copyToDeviceTime = timer.elapsed();
    std::cout << "Initialize rand state and copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';
#endif
//Main time loop
    for(int tt=0; tt< nT; tt++)
    {
#ifdef __CUDACC__
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), 
                move_boris(dt,BoundaryDevicePointer, nLines,
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer,
                    nR_PreSheathEfield,nZ_PreSheathEfield,EfieldGridRDevicePointer,
                    EfieldGridZDevicePointer,EfieldRDevicePointer,
                    EfieldZDevicePointer,EfieldTDevicePointer));
        
        
        try {
            thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    geometry_check(nLines,BoundaryDevicePointer,dt,tt) );
        }
        catch (thrust::system_error &e) {
            std::cerr << "Thrust system error: " << e.what() << std::endl;
            exit(-1);
        }
#if USEIONIZATION > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                ionize(dt,
                    nR_Dens,nZ_Dens,DensGridRDevicePointer,DensGridZDevicePointer ,NeDevicePointer,  
                    nR_Temp,nZ_Temp,TempGridRDevicePointer,TempGridZDevicePointer,TeDevicePointer,
                    nTemperaturesIonize, nDensitiesIonize,IonizationTemperatureGridDevicePointer, 
                    IonizationDensityGridDevicePointer,IonizationRateCoeffDevicePointer));
#endif
#if USERECOMBINATION > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                crossFieldDiffusion(dt,perpDiffusionCoeff,
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));
            
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                    geometry_check(nLines,BoundaryDevicePointer,dt,tt) );
#endif
#if USECOULOMBCOLLISIONS > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), 
                coulombCollisions(dt,
                    nR_flowV,nZ_flowV,FlowVGridRDevicePointer,FlowVGridZDevicePointer,
                    FlowVrDevicePointer,FlowVzDevicePointer,FlowVtDevicePointer,
                    nR_Dens,nZ_Dens,DensGridRDevicePointer,DensGridZDevicePointer,
                    NeDevicePointer,    
                    nR_Temp,nZ_Temp,TempGridRDevicePointer,TempGridZDevicePointer,
                    TeDevicePointer,
                    background_Z,background_amu, 
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));

#endif
#if USETHERMALFORCE > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(),
                thermalForce(dt,background_amu,
                    nR_gradT,nZ_gradT,GradTGridRDevicePointer,GradTGridZDevicePointer,
                    GradTiRDevicePointer,GradTiZDevicePointer, GradTiTDevicePointer, 
                    GradTeRDevicePointer, GradTeZDevicePointer, GradTeTDevicePointer, 
                    nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
                    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));
#endif

#if USESURFACEMODEL > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), 
                reflection(dt,nLines,BoundaryDevicePointer) );
#endif        
#else
#if USE_BOOST
cpu_times moveTime0 = timer.elapsed();
#endif
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
                move_boris(dt,hostBoundaryVector.data(),nLines, 
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
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
            geometry_check(nLines,hostBoundaryVector.data(),dt,tt) );
#if USE_BOOST
cpu_times geomTime1 = timer.elapsed();
geomCheckTime = geomCheckTime + (geomTime1.wall - geomTime0.wall);
#endif
#if USEIONIZATION > 0
#if USE_BOOST
cpu_times ionizTime0 = timer.elapsed();
#endif
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(),
            ionize(dt,
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
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
            recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
                crossFieldDiffusion(dt,perpDiffusionCoeff,
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front()));
        
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
                geometry_check(nLines,hostBoundaryVector.data(),dt,tt) );
#endif
#if USECOULOMBCOLLISIONS > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
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
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
                thermalForce(dt,background_amu,
                    nR_gradT,nZ_gradT,&gradTGridr.front(),&gradTGridz.front(),
                    &gradTiR.front(),&gradTiZ.front(),&gradTiT.front(),
                    &gradTeR.front(),&gradTeZ.front(),&gradTeT.front(),
                    nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                    &br.front(),&bz.front(),&bt.front()));
#endif
#if USESURFACEMODEL > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), 
                reflection(dt,nLines,hostBoundaryVector.data()) );
#endif        
#if PARTICLE_TRACKS >0
if (tt % subSampleFac == 0)  
{    
        for(int i=0;i<nP;i++)
        {
            positionHistoryX[i][tt/subSampleFac] = hostCudaParticleVector[i].xprevious;
            positionHistoryY[i][tt/subSampleFac] = hostCudaParticleVector[i].yprevious;
            positionHistoryZ[i][tt/subSampleFac] = hostCudaParticleVector[i].zprevious;
            velocityHistoryX[i][tt/subSampleFac] = hostCudaParticleVector[i].vx;
            velocityHistoryY[i][tt/subSampleFac] = hostCudaParticleVector[i].vy;
            velocityHistoryZ[i][tt/subSampleFac] = hostCudaParticleVector[i].vz;
        }
}
#endif
#endif
    }
#if USE_BOOST
    cpu_times ionizeTimeGPU = timer.elapsed();
    std::cout << "Particle Moving Time: " << ionizeTimeGPU.wall*1e-9 << '\n';
#endif
#ifdef __CUDACC__
    hostCudaParticleVector = deviceCudaParticleVector;
#endif

int ring1 = 0;
int ring2 = 0;
int noWall = 0;
float meanTransitTime = 0.0;

for(int i=0; i<nP ; i++)
{
	if(hostCudaParticleVector[i].wallIndex == boundaryIndex_ImpurityLaunch[0])
	{
		ring1++;
	}
	else if(hostCudaParticleVector[i].wallIndex == boundaryIndex_ImpurityLaunch[1])
	{
		ring2++;
	}
	
	if(hostCudaParticleVector[i].wallIndex == 0)
	{
		noWall++;
	}
	
	meanTransitTime = meanTransitTime + hostCudaParticleVector[i].transitTime;
	
} 
meanTransitTime = meanTransitTime/(nP-noWall);
std::cout << "Number of impurity particles deposited on ring 1 " << ring1 << std::endl;
std::cout << "Number of impurity particles deposited on ring 2 " << ring2 << std::endl;
std::cout << "Number of impurity particles not deposited " << noWall << std::endl;
std::cout << "Mean transit time of deposited particles " << meanTransitTime << std::endl;
   //for(int i=0; i < hostCudaParticleVector.size(); i++){
       //std::cout << " final pos" <<  i << " " <<hostCudaParticleVector[i].x << " " << hostCudaParticleVector[i].y << " " << hostCudaParticleVector[i].z << std::endl;
        /*if(hostCudaParticleVector[i].hitWall == 1){
        surfaceIndexY = int(floor((hostCudaParticleVector[i].y - yMin)/(yMax - yMin)*(nY) + 0.0f));
        surfaceIndexZ = int(floor((hostCudaParticleVector[i].z - zMin)/(zMax - zMin)*(nZ) + 0.0f));
        SurfaceBins[surfaceIndexY][surfaceIndexZ] +=  1.0 ;

        SurfaceBinsCharge[surfaceIndexY][surfaceIndexZ] += hostCudaParticleVector[i].Z ;
        SurfaceBinsEnergy[surfaceIndexY][surfaceIndexZ] += 0.5*hostCudaParticleVector[i].amu*1.6737236e-27*(hostCudaParticleVector[i].vx*hostCudaParticleVector[i].vx +  hostCudaParticleVector[i].vy*hostCudaParticleVector[i].vy+ hostCudaParticleVector[i].vz*hostCudaParticleVector[i].vz)/1.60217662e-19;
        } */ 
   // }

//    OUTPUT( outname,nY, nZ, SurfaceBins);
//    OUTPUT( outnameCharge,nY, nZ, SurfaceBinsCharge);
//    OUTPUT( outnameEnergy,nY, nZ, SurfaceBinsEnergy);

    ofstream outfile2;
    outfile2.open ("positions.m");
    for(int i=1 ; i<=nP ; i++)
      {
        outfile2 << "Pos( " << i<< ",:) = [ " ;
        outfile2 << hostCudaParticleVector[i-1].x << " " << hostCudaParticleVector[i-1].y 
            << " " << hostCudaParticleVector[i-1].z << " ];" << std::endl;
      }
       outfile2.close();
// Write netCDF output for positions
for (int i=0; i<nP; i++)
{
    finalPosX[i] = hostCudaParticleVector[i].xprevious;
    finalPosY[i] = hostCudaParticleVector[i].yprevious;
    finalPosZ[i] = hostCudaParticleVector[i].zprevious;
    finalVx[i] = hostCudaParticleVector[i].vx;
    finalVy[i] = hostCudaParticleVector[i].vy;
    finalVz[i] = hostCudaParticleVector[i].vz;
    transitTime[i] = hostCudaParticleVector[i].transitTime;
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
#if PARTICLE_TRACKS > 0
/*char outnameX[] = "positionHistoryX.m";
OUTPUT( outnameX,nP, nT/subSampleFac, positionHistoryX);
char outnameY[] = "positionHistoryY.m";
OUTPUT( outnameY,nP, nT/subSampleFac, positionHistoryY);
char outnameZ[] = "positionHistoryZ.m";
OUTPUT( outnameZ,nP, nT/subSampleFac, positionHistoryZ);
char outnameVX[] = "velocityHistoryX.m";
OUTPUT( outnameVX,nP, nT/subSampleFac,velocityHistoryX);
char outnameVY[] = "velocityHistoryY.m";
OUTPUT( outnameVY,nP, nT/subSampleFac, velocityHistoryY);
char outnameVZ[] = "velocityHistoryZ.m";
OUTPUT( outnameVZ,nP, nT/subSampleFac, velocityHistoryZ);
*/

// Write netCDF output for histories
NcFile ncFile("history.nc", NcFile::replace);
NcDim nc_nT = ncFile.addDim("nT",nT/subSampleFac);
NcDim nc_nP = ncFile.addDim("nP",nP);
vector<NcDim> dims;
dims.push_back(nc_nP);
dims.push_back(nc_nT);

NcVar nc_x = ncFile.addVar("x",ncDouble,dims);
NcVar nc_y = ncFile.addVar("y",ncDouble,dims);
NcVar nc_z = ncFile.addVar("z",ncDouble,dims);

NcVar nc_vx = ncFile.addVar("vx",ncDouble,dims);
NcVar nc_vy = ncFile.addVar("vy",ncDouble,dims);
NcVar nc_vz = ncFile.addVar("vz",ncDouble,dims);

nc_x.putVar(positionHistoryX[0]);
nc_y.putVar(positionHistoryY[0]);
nc_z.putVar(positionHistoryZ[0]);

nc_vx.putVar(velocityHistoryX[0]);
nc_vy.putVar(velocityHistoryY[0]);
nc_vz.putVar(velocityHistoryZ[0]);

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
    return 0;
}
