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
#include <boost/timer/timer.hpp>
#include "boost/filesystem.hpp"
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
using namespace boost::timer;
using namespace netCDF;
using namespace exceptions;

int main()
{

Config cfg,cfg_geom;

cfg.readFile("gitrInput.cfg");
cfg_geom.readFile("gitrGeometry.cfg");


// Background species info
double background_Z = cfg.lookup("backgroundPlasmaProfiles.Z");
double background_amu = cfg.lookup("backgroundPlasmaProfiles.amu");

std::string outnameBfieldR = "BfieldR.m";
std::string outnameBfieldZ = "BfieldZ.m";
std::string outnameBfieldT = "BfieldT.m";
std::string outnameGridR = "gridR.m";
std::string outnameGridZ = "gridZ.m";
std::string profiles_folder = "profiles";
#if BFIELD_INTERP == 0
int nR_Bfield = 1;
int nZ_Bfield = 1;
std::vector<double> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
std::vector<double> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);
br[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.br");
bz[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bz");
bt[0] = cfg.lookup("backgroundPlasmaProfiles.Bfield.bt");
#else
int nR_Bfield;
int nZ_Bfield;

int b1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Bfield.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNrString"),
            cfg.lookup("backgroundPlasmaProfiles.Bfield.gridNzString"),nR_Bfield,nZ_Bfield);

std::vector<double> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
std::vector<double> br(nR_Bfield*nZ_Bfield), bz(nR_Bfield*nZ_Bfield),bt(nR_Bfield*nZ_Bfield);

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

OUTPUT1d(profiles_folder,outnameGridR, nR_Bfield, &bfieldGridr.front());
OUTPUT1d(profiles_folder,outnameGridZ, nZ_Bfield, &bfieldGridz.front());
OUTPUT2d(profiles_folder,outnameBfieldR, nR_Bfield, nZ_Bfield, &br.front());
OUTPUT2d(profiles_folder,outnameBfieldZ, nR_Bfield, nZ_Bfield, &bz.front());
OUTPUT2d(profiles_folder,outnameBfieldT, nR_Bfield, nZ_Bfield, &bt.front());
#endif
#ifdef __CUDACC__
    thrust::device_vector<double> deviceBfieldGridRVector = bfieldGridr;
    thrust::device_vector<double> deviceBfieldGridZVector = bfieldGridz;
    thrust::device_vector<double> deviceBfieldRVector = br;
    thrust::device_vector<double> deviceBfieldZVector = bz;
    thrust::device_vector<double> deviceBfieldTVector = bt;

    double * BfieldGridRDevicePointer = thrust::raw_pointer_cast(deviceBfieldGridRVector.data());
    double * BfieldGridZDevicePointer = thrust::raw_pointer_cast(deviceBfieldGridZVector.data());
    double * BfieldRDevicePointer = thrust::raw_pointer_cast(deviceBfieldRVector.data());
    double * BfieldZDevicePointer = thrust::raw_pointer_cast(deviceBfieldZVector.data());
    double * BfieldTDevicePointer = thrust::raw_pointer_cast(deviceBfieldTVector.data());
#endif
#if TEMP_INTERP == 0
int nR_Temp = 1;
int nZ_Temp = 1;
std::vector<double> TempGridr(nR_Temp), TempGridz(nZ_Temp);
std::vector<double> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);
ti[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.ti");
te[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.te");
#else
int nR_Temp;
int nZ_Temp;

int t1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNrString"),
            cfg.lookup("backgroundPlasmaProfiles.Temperature.gridNzString"),nR_Temp,nZ_Temp);
std::cout << "t1 " << nR_Temp << nZ_Temp << std::endl;
std::vector<double> TempGridr(nR_Temp), TempGridz(nZ_Temp);
std::vector<double> ti(nR_Temp*nZ_Temp), te(nR_Temp*nZ_Temp);

int t2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Temperature.gridRString"), TempGridr);

int t3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Temperature.gridZString"), TempGridz);
std::cout << "temperature import" << nZ_Temp << nR_Temp << std::endl;
int t4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Temperature.IonTempString"), ti);

int t5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Temperature.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Temperature.ElectronTempString"), te);
std::string outnameTi = "ti.m";
std::string outnameTe = "te.m";
OUTPUT2d(profiles_folder,outnameTi, nR_Temp, nZ_Temp, &ti.front());
OUTPUT2d(profiles_folder,outnameTe, nR_Temp, nZ_Temp, &te.front());
#endif

#ifdef __CUDACC__
    thrust::device_vector<double> deviceTempGridRVector = TempGridr;
    thrust::device_vector<double> deviceTempGridZVector = TempGridz;
    thrust::device_vector<double> deviceTiVector = ti;
    thrust::device_vector<double> deviceTeVector = te;
    double * TempGridRDevicePointer = thrust::raw_pointer_cast(deviceTempGridRVector.data());
    double * TempGridZDevicePointer = thrust::raw_pointer_cast(deviceTempGridZVector.data());
    double * TiDevicePointer = thrust::raw_pointer_cast(deviceTiVector.data());
    double * TeDevicePointer = thrust::raw_pointer_cast(deviceTeVector.data());
#endif


#if DENSITY_INTERP == 0
int nR_Dens = 1;
int nZ_Dens = 1;
std::vector<double> DensGridr(nR_Dens), DensGridz(nZ_Dens);
std::vector<double> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
ni[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ni");
ne[0] = cfg.lookup("backgroundPlasmaProfiles.Density.ne");
#else
int nR_Dens;
int nZ_Dens;

int n1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.gridNrString"),
            cfg.lookup("backgroundPlasmaProfiles.Density.gridNzString"),nR_Dens,nZ_Dens);

std::vector<double> DensGridr(nR_Dens), DensGridz(nZ_Dens);
std::vector<double> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);

int n2 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Density.gridRString"), DensGridr);

int n3 = read_profile1d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Density.gridZString"), DensGridz);

int n4 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Density.IonDensityString"), ni);

int n5 = read_profile2d(cfg.lookup("backgroundPlasmaProfiles.Density.fileString"),
        cfg.lookup("backgroundPlasmaProfiles.Density.ElectronDensityString"), ne);
std::string outnameNi = "ni.m";
std::string outnameNe = "ne.m";
OUTPUT2d(profiles_folder,outnameNi, nR_Dens, nZ_Dens, &ni.front());
OUTPUT2d(profiles_folder,outnameNe, nR_Dens, nZ_Dens, &ne.front());
#endif

#ifdef __CUDACC__
    thrust::device_vector<double> deviceDensGridRVector = DensGridr;
    thrust::device_vector<double> deviceDensGridZVector = DensGridz;
    thrust::device_vector<double> deviceNiVector = ni;
    thrust::device_vector<double> deviceNeVector = ne;
    double * DensGridRDevicePointer = thrust::raw_pointer_cast(deviceDensGridRVector.data());
    double * DensGridZDevicePointer = thrust::raw_pointer_cast(deviceDensGridZVector.data());
    double * NiDevicePointer = thrust::raw_pointer_cast(deviceNiVector.data());
    double * NeDevicePointer = thrust::raw_pointer_cast(deviceNeVector.data());
#endif

#if FLOWV_INTERP == 0
int nR_flowV = 1;
int nZ_flowV = 1;
std::vector<double> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
std::vector<double> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),flowVt(nR_flowV*nZ_flowV);
flowVr[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVr");
flowVz[0] = cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.flowVz");
#else
int nR_flowV;
int nZ_flowV;

int f1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNrString"),
            cfg.lookup("backgroundPlasmaProfiles.FlowVelocity.gridNzString"),nR_flowV,nZ_flowV);

std::vector<double> flowVGridr(nR_flowV), flowVGridz(nZ_flowV);
std::vector<double> flowVr(nR_flowV*nZ_flowV), flowVz(nR_flowV*nZ_flowV),flowVt(nR_flowV*nZ_flowV);

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

std::string outnameFlowVr = "flowVr.m";
std::string outnameFlowVz = "flowVz.m";
std::string outnameFlowVt = "flowVt.m";
OUTPUT2d(profiles_folder,outnameFlowVr, nR_flowV, nZ_flowV, &flowVr.front());
OUTPUT2d(profiles_folder,outnameFlowVz, nR_flowV, nZ_flowV, &flowVz.front());
OUTPUT2d(profiles_folder,outnameFlowVt, nR_flowV, nZ_flowV, &flowVt.front());
#endif

#ifdef __CUDACC__
    thrust::device_vector<double> deviceFlowVGridRVector = flowVGridr;
    thrust::device_vector<double> deviceFlowVGridZVector = flowVGridz;
    thrust::device_vector<double> deviceFlowVrVector = flowVr;
    thrust::device_vector<double> deviceFlowVzVector = flowVz;
    thrust::device_vector<double> deviceFlowVtVector = flowVt;
    double * FlowVGridRDevicePointer = thrust::raw_pointer_cast(deviceFlowVGridRVector.data());
    double * FlowVGridZDevicePointer = thrust::raw_pointer_cast(deviceFlowVGridZVector.data());
    double * FlowVrDevicePointer = thrust::raw_pointer_cast(deviceFlowVrVector.data());
    double * FlowVzDevicePointer = thrust::raw_pointer_cast(deviceFlowVzVector.data());
    double * FlowVtDevicePointer = thrust::raw_pointer_cast(deviceFlowVtVector.data());
#endif

#if GRADT_INTERP == 0
int nR_gradT = 1;
int nZ_gradT = 1;
std::vector<double> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
std::vector<double> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),
    gradTiR(nR_gradT*nZ_gradT), gradTiZ(nR_gradT*nZ_gradT);
gradTeR[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeR");
gradTeZ[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTeZ");
gradTiR[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiR");
gradTiZ[0] = cfg.lookup("backgroundPlasmaProfiles.gradT.gradTiZ");
#else
int nR_gradT;
int nZ_gradT;

int g1 = read_profileNs(cfg.lookup("backgroundPlasmaProfiles.gradT.fileString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gridNrString"),
            cfg.lookup("backgroundPlasmaProfiles.gradT.gridNzString"),nR_gradT,nZ_gradT);

std::vector<double> gradTGridr(nR_gradT), gradTGridz(nZ_gradT);
std::vector<double> gradTeR(nR_gradT*nZ_gradT), gradTeZ(nR_gradT*nZ_gradT),
    gradTiR(nR_gradT*nZ_gradT), gradTiZ(nR_gradT*nZ_gradT);

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
std::string outnameGradTiR = "gradTiR.m";
std::string outnameGradTiZ = "gradTiZ.m";
std::string outnameGradTeR = "gradTeR.m";
std::string outnameGradTeZ = "gradTeZ.m";
OUTPUT2d(profiles_folder,outnameGradTiR, nR_gradT, nZ_gradT, &gradTiR.front());
OUTPUT2d(profiles_folder,outnameGradTiZ, nR_gradT, nZ_gradT, &gradTiZ.front());
OUTPUT2d(profiles_folder,outnameGradTeR, nR_gradT, nZ_gradT, &gradTeR.front());
OUTPUT2d(profiles_folder,outnameGradTeZ, nR_gradT, nZ_gradT, &gradTeZ.front());
#endif

#ifdef __CUDACC__
    thrust::device_vector<double> deviceGradTGridRVector = gradTGridr;
    thrust::device_vector<double> deviceGradTGridZVector = gradTGridz;
    thrust::device_vector<double> deviceGradTiRVector = gradTiR;
    thrust::device_vector<double> deviceGradTiZVector = gradTiZ;
    thrust::device_vector<double> deviceGradTeRVector = gradTeR;
    thrust::device_vector<double> deviceGradTeZVector = gradTeZ;
    double * GradTGridRDevicePointer = thrust::raw_pointer_cast(deviceGradTGridRVector.data());
    double * GradTGridZDevicePointer = thrust::raw_pointer_cast(deviceGradTGridZVector.data());
    double * GradTiRDevicePointer = thrust::raw_pointer_cast(deviceGradTiRVector.data());
    double * GradTiZDevicePointer = thrust::raw_pointer_cast(deviceGradTiZVector.data());
    double * GradTeRDevicePointer = thrust::raw_pointer_cast(deviceGradTeRVector.data());
    double * GradTeZDevicePointer = thrust::raw_pointer_cast(deviceGradTeZVector.data());
#endif

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

std::vector<double> rateCoeff_Ionization(nCS_Ionize*nTemperaturesIonize*nDensitiesIonize);
std::vector<double> gridTemperature_Ionization(nTemperaturesIonize),
    gridDensity_Ionization(nDensitiesIonize);

int    i2 = read_profiles(cfg.lookup("impurityParticleSource.ionization.fileString"),
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

std::vector<double> rateCoeff_Recombination(nCS_Recombine*nTemperaturesRecombine*nDensitiesRecombine);
std::vector<double> gridTemperature_Recombination(nTemperaturesRecombine),
    gridDensity_Recombination(nDensitiesRecombine);

int    i4 = read_profiles(cfg.lookup("impurityParticleSource.recombination.fileString"),
        nTemperaturesRecombine,nDensitiesRecombine,
        cfg.lookup("impurityParticleSource.recombination.TempGridVarName"), 
        gridTemperature_Recombination,cfg.lookup("impurityParticleSource.recombination.DensGridVarName"),
        gridDensity_Recombination,
        cfg.lookup("impurityParticleSource.recombination.CoeffVarName"),
        rateCoeff_Recombination);

#ifdef __CUDACC__
    thrust::device_vector<double> deviceTemperatureGridIonizationVector = gridTemperature_Ionization;
    thrust::device_vector<double> deviceDensityGridIonizationVector = gridDensity_Ionization;
    thrust::device_vector<double> deviceTemperatureGridRecombinationVector = gridTemperature_Recombination;
    thrust::device_vector<double> deviceDensityGridRecombinationVector = gridDensity_Recombination;
    thrust::device_vector<double> deviceRateCoeffIonizationVector = rateCoeff_Ionization;
    thrust::device_vector<double> deviceRateCoeffRecombinationVector = rateCoeff_Recombination;
    double * IonizationTemperatureGridDevicePointer = thrust::raw_pointer_cast(deviceTemperatureGridIonizationVector.data());
    double * RecombinationTemperatureGridDevicePointer = thrust::raw_pointer_cast(deviceTemperatureGridRecombinationVector.data());
    double * IonizationDensityGridDevicePointer = thrust::raw_pointer_cast(deviceDensityGridIonizationVector.data());
    double * RecombinationDensityGridDevicePointer = thrust::raw_pointer_cast(deviceDensityGridRecombinationVector.data());
    double * IonizationRateCoeffDevicePointer = thrust::raw_pointer_cast(deviceRateCoeffIonizationVector.data());
    double * RecombinationRateCoeffDevicePointer = thrust::raw_pointer_cast(deviceRateCoeffRecombinationVector.data());
#endif

char outname[] = "Deposition.m";
char outnameCharge[] = "Charge.m";
char outnameEnergy[] = "Energy.m";

//Geometry Definition
Setting& geom = cfg_geom.lookup("geom");
int nLines = geom["x1"].getLength();
std::cout << "Number of Geometric Objects Loaded: " << nLines << std::endl;

std::vector<Boundary> hostBoundaryVector(nLines+1);

std::string geom_outname = "geom.m";
std::string geom_folder = "geometry";
       ofstream outfile;
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
       std::string full_path = geom_folder + "/" + geom_outname;
      outfile.open (full_path );

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
std::for_each(hostBoundaryVector.begin(), hostBoundaryVector.end()-1, boundary_init(background_Z,background_amu,nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ni.front(),nR_Bfield,nZ_Bfield,&bfieldGridr.front(),&bfieldGridz.front(),&br.front(),&bz.front(), &bt.front(),
       nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&ti.front() ));

#ifdef __CUDACC__
    thrust::device_vector<Boundary> deviceBoundaryVector = hostBoundaryVector;
    Boundary * BoundaryDevicePointer = thrust::raw_pointer_cast(deviceBoundaryVector.data());
#else
    std::vector<Boundary> * BoundaryHostPointer = &hostBoundaryVector;    
#endif
//Efield
    
std::vector<double> Efieldr(nR_Bfield*nZ_Bfield), Efieldz(nR_Bfield*nZ_Bfield),Efieldt(nR_Bfield*nZ_Bfield),minDist(nR_Bfield*nZ_Bfield);
double thisE[3] = {0.0,0.0,0.0};

for(int i=0;i<nR_Bfield;i++)
{
    for(int j=0;j<nZ_Bfield;j++)
    {
        minDist[(nR_Bfield - 1 -i)*nZ_Bfield+(nZ_Bfield -1-j)] = getE ( bfieldGridr[i], 0.0, bfieldGridz[j],  thisE, hostBoundaryVector,nLines );
        //std::cout << "point " << bfieldGridr[i] << " " << bfieldGridz[j] << " dist" << minDist[(i)*nZ_Bfield+(j)] << std::endl;
        Efieldr[i*nZ_Bfield+j] = thisE[0];
        Efieldz[i*nZ_Bfield+j] = thisE[2];
        Efieldt[i*nZ_Bfield+j] = thisE[1];
    }
}
for(int i=0;i<nR_Bfield;i++)
{
    for(int j=0;j<nZ_Bfield;j++)
    {
       // std::cout << "point " << bfieldGridr[i] << " " << bfieldGridz[j] << " dist" << minDist[(i)+(j)*nR_Bfield] << std::endl;
    }
}
std::string outnameEfieldR = "EfieldR.m";
std::string outnameEfieldZ = "EfieldZ.m";
std::string outnameEfieldT = "EfieldT.m";
std::string outnameMinDist = "DistToSurface.m";
OUTPUT2d(profiles_folder,outnameEfieldR, nR_Bfield, nZ_Bfield, &Efieldr.front());
OUTPUT2d(profiles_folder,outnameEfieldZ, nR_Bfield, nZ_Bfield, &Efieldz.front());
OUTPUT2d(profiles_folder,outnameEfieldT, nR_Bfield, nZ_Bfield, &Efieldt.front());
OUTPUT2d(profiles_folder,outnameMinDist, nR_Bfield, nZ_Bfield, &minDist.front());

// Particle time stepping control

int ionization_nDtPerApply  = cfg.lookup("timeStep.ionization_nDtPerApply");
int collision_nDtPerApply  = cfg.lookup("timeStep.collision_nDtPerApply");
// Perp DiffusionCoeff - only used when Diffusion interpolator is = 0
double perpDiffusionCoeff = cfg.lookup("backgroundPlasmaProfiles.Diffusion.Dperp");

// Background species info
//double background_Z = cfg.lookup("backgroundPlasmaProfiles.Z");
//double background_amu = cfg.lookup("backgroundPlasmaProfiles.amu");
double *background_flow;
double *maxDensity;
double *maxTemp_eV;

#ifdef __CUDACC__
    cout<<"Using THRUST"<<endl;
#else
    cout<<"Not using THRUST"<<endl;
#endif

/*Setting& backgroundPlasma = cfg.lookup("backgroundPlasma");
int nS = backgroundPlasma["Z"].getLength();

Setting& diagnostics = cfg.lookup("diagnostics");


background_Z = new int[nS];
background_amu = new double[nS];
background_flow = new double[nS];
maxDensity = new double[nS];
maxTemp_eV = new double[nS];

for(int i=0; i<nS; i++)
{
background_Z[i] = backgroundPlasma["Z"][i];
background_amu[i] = backgroundPlasma["amu"][i];
background_flow[i] = backgroundPlasma["flow"]["fractionOfThermalVelocity"][i];
maxDensity[i] = backgroundPlasma["density"]["max"][i];
maxTemp_eV[i] = backgroundPlasma["temp"]["max"][i];
}
*/
    double x = cfg.lookup("impurityParticleSource.initialConditions.x_start");
    double y = cfg.lookup("impurityParticleSource.initialConditions.y_start");
    double z = cfg.lookup("impurityParticleSource.initialConditions.z_start");
    
    double Ex = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
    double Ey = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
    double Ez = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");
    
    double amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
    double Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");
    double charge = cfg.lookup("impurityParticleSource.initialConditions.charge");
    /*
    double **SurfaceBins;
    double **SurfaceBinsCharge;
    double **SurfaceBinsEnergy;
    double **SurfaceBinsErosion;
    
    SurfaceBins = new double*[nY];
    SurfaceBinsCharge = new double*[nY];
    SurfaceBinsEnergy = new double*[nY];
    SurfaceBinsErosion = new double*[nY];

    SurfaceBins[0] = new double[nY*nZ];
    SurfaceBinsCharge[0] = new double[nY*nZ];
    SurfaceBinsEnergy[0] = new double[nY*nZ];
    SurfaceBinsErosion[0] = new double[nY*nZ];
            
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
    double dt;
    double nPtsPerGyroOrbit = cfg.lookup("timeStep.nPtsPerGyroOrbit");
    dt = 1e-6/nPtsPerGyroOrbit;

    int nP = cfg.lookup("impurityParticleSource.nP");
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
    double impurity_Z = cfg.lookup("impurityParticleSource.Z");
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
    double randDevice = 6.5298E+5;
    std::default_random_engine generator0(randDevice);
#endif
    
    std::vector<float> boundarySeeds0(4*nImpurityBoundaries);
    std::generate( boundarySeeds0.begin(), boundarySeeds0.end(), [&]() { return distributionForSeeds(generator0); } );
    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    double rand0 = 0.0;
    double rand1 = 0.0;
    double rand2 = 0.0;
    double rand3 = 0.0;

    std::vector<std::mt19937> s0(4*nImpurityBoundaries);
    
    double E0 = 0.0;
//Create Thompson Distribution
    double surfaceBindingEnergy = cfg.lookup("impurityParticleSource.source_material_SurfaceBindingEnergy");
    std::cout << "surface binding energy " << surfaceBindingEnergy << std::endl;
    int nThompDistPoints = 200;
    double max_Energy = 100.0;
    std::vector<double> ThompsonDist(nThompDistPoints),CumulativeDFThompson(nThompDistPoints);
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
            x = hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].x1 + hostBoundaryVector[boundaryIndex_ImpurityLaunch[i]].length*rand0;//1.4290;
            //std::cout << "start pos 1 " << x << std::endl;
            z = -1.2540+0.0001;
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
            z = -1.3660+0.0001;
            rand1 = dist01(s0[5]);
            rand2 = dist01(s0[6]);
            rand3 = dist01(s0[7]);
            E0 = interp1dUnstructured(rand2,nThompDistPoints, max_Energy, &CumulativeDFThompson.front());
            Ex = E0*cos(3.1415*rand1)*sin(3.1415*rand3);
            Ey = E0*cos(3.1415*rand3);
            Ez = E0*sin(3.1415*rand1)*sin(3.1415*rand3);
        }
            
        Particle p1(x,0.0,z,Ex,0.0,Ez,74,184.0,charge);
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
            {   double theta = dist2(rd2)*2*3.1415;
                double phi = dist2(rd2)*3.1415;
                double mag = 2e3;
                hostCudaParticleVector[i].vx = mag*cos(theta)*sin(phi);
                hostCudaParticleVector[i].vy = mag*sin(theta)*sin(phi);
                hostCudaParticleVector[i].vz = mag*cos(phi);
            }
#endif
#if PARTICLE_TRACKS > 0
int subSampleFac = 10;
double **positionHistoryX;
double **positionHistoryY;
double **positionHistoryZ;
double **velocityHistoryX;
double **velocityHistoryY;
double **velocityHistoryZ;
positionHistoryX = new double* [nP];
positionHistoryY = new double* [nP];
positionHistoryZ = new double* [nP];
velocityHistoryX = new double* [nP];
velocityHistoryY = new double* [nP];
velocityHistoryZ = new double* [nP];
positionHistoryX[0] = new double [nT*nP/subSampleFac];
positionHistoryY[0] = new double [nT*nP/subSampleFac];
positionHistoryZ[0] = new double [nT*nP/subSampleFac];
velocityHistoryX[0] = new double [nT*nP/subSampleFac];
velocityHistoryY[0] = new double [nT*nP/subSampleFac];
velocityHistoryZ[0] = new double [nT*nP/subSampleFac];
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
            cpu_timer timer;

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
    double ionization_seeds = cfg.lookup("operators.ionization.seed");
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
    double moveTime = 0.0;
    double geomCheckTime = 0.0;
    double ionizTime = 0.0;
    cpu_times copyToDeviceTime = timer.elapsed();
    std::cout << "Initialize rand state and copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';
    for(int tt=0; tt< nT; tt++)
    {
#ifdef __CUDACC__
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), move_boris(dt,BoundaryDevicePointer, nLines,nR_Bfield,nZ_Bfield, BfieldGridRDevicePointer,BfieldGridZDevicePointer,
    BfieldRDevicePointer,BfieldZDevicePointer,BfieldTDevicePointer));
        try {
            thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), geometry_check(nLines,BoundaryDevicePointer) );
        }
        catch (thrust::system_error &e) {
            std::cerr << "Thrust system error: " << e.what() << std::endl;
            exit(-1);
        }
#if USEIONIZATION > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), ionize(dt) );
#endif
#if USERECOMBINATION > 0
    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), crossFieldDiffusion(dt,perpDiffusionCoeff));
#endif
#if USECOULOMBCOLLISIONS > 0
    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), coulombCollisions(dt) );
#endif
#if USETHERMALFORCE > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), thermalForce(dt) );
#endif
#else
cpu_times moveTime0 = timer.elapsed();
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), move_boris(dt,hostBoundaryVector,nLines, nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                &br.front(),&bz.front(),&bt.front()));
cpu_times moveTime1 = timer.elapsed();
moveTime = moveTime + (moveTime1.wall - moveTime0.wall);

cpu_times geomTime0 = timer.elapsed();
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), geometry_check(nLines,hostBoundaryVector) );
cpu_times geomTime1 = timer.elapsed();
geomCheckTime = geomCheckTime + (geomTime1.wall - geomTime0.wall);
#if USEIONIZATION > 0
cpu_times ionizTime0 = timer.elapsed();
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), ionize(dt,
                nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),
                nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                nTemperaturesIonize, nDensitiesIonize, &gridTemperature_Ionization.front(),
               &gridDensity_Ionization.front(), &rateCoeff_Ionization.front() ) );
cpu_times ionizTime1 = timer.elapsed();
ionizTime = ionizTime + (ionizTime1.wall - ionizTime0.wall);
#endif
#if USERECOMBINATION > 0
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
    //std::cout<< "Perp diffusion loop " << perpDiffusionCoeff << std::endl;
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), crossFieldDiffusion(dt,perpDiffusionCoeff,nR_Bfield,nZ_Bfield, &bfieldGridr.front(),&bfieldGridz.front(),
                                    &br.front(),&bz.front(),&bt.front()));
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), geometry_check(nLines,hostBoundaryVector) );
#endif
#if USECOULOMBCOLLISIONS > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), coulombCollisions(dt,nR_flowV,nZ_flowV,&flowVGridr.front(),&flowVGridz.front(),&flowVr.front(),&flowVz.front(),
                    &flowVt.front(),
                nR_Dens,nZ_Dens,&DensGridr.front(),&DensGridz.front(),&ne.front(),
                nR_Temp,nZ_Temp,&TempGridr.front(),&TempGridz.front(),&te.front(),
                background_Z,background_amu,nR_Bfield,nZ_Bfield, &bfieldGridr.front(),
                &bfieldGridz.front(),&br.front(),&bz.front(),&bt.front()));
#endif
#if USETHERMALFORCE > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), thermalForce(dt,background_amu,nR_gradT,nZ_gradT,&gradTGridr.front(),&gradTGridz.front(),
                    &gradTiR.front(),&gradTiZ.front(),&gradTeR.front(),&gradTeZ.front() ) );
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
    cpu_times ionizeTimeGPU = timer.elapsed();
    std::cout << "Particle Moving Time: " << ionizeTimeGPU.wall*1e-9 << '\n';

#ifdef __CUDACC__
    hostCudaParticleVector = deviceCudaParticleVector;
#endif

    for(int i=0; i < hostCudaParticleVector.size(); i++){
       //std::cout << " final pos" <<  i << " " <<hostCudaParticleVector[i].x << " " << hostCudaParticleVector[i].y << " " << hostCudaParticleVector[i].z << std::endl;
        /*if(hostCudaParticleVector[i].hitWall == 1){
        surfaceIndexY = int(floor((hostCudaParticleVector[i].y - yMin)/(yMax - yMin)*(nY) + 0.0f));
        surfaceIndexZ = int(floor((hostCudaParticleVector[i].z - zMin)/(zMax - zMin)*(nZ) + 0.0f));
        SurfaceBins[surfaceIndexY][surfaceIndexZ] +=  1.0 ;

        SurfaceBinsCharge[surfaceIndexY][surfaceIndexZ] += hostCudaParticleVector[i].Z ;
        SurfaceBinsEnergy[surfaceIndexY][surfaceIndexZ] += 0.5*hostCudaParticleVector[i].amu*1.6737236e-27*(hostCudaParticleVector[i].vx*hostCudaParticleVector[i].vx +  hostCudaParticleVector[i].vy*hostCudaParticleVector[i].vy+ hostCudaParticleVector[i].vz*hostCudaParticleVector[i].vz)/1.60217662e-19;
        } */ 
    }

//    OUTPUT( outname,nY, nZ, SurfaceBins);
//    OUTPUT( outnameCharge,nY, nZ, SurfaceBinsCharge);
//    OUTPUT( outnameEnergy,nY, nZ, SurfaceBinsEnergy);

    ofstream outfile2;
    outfile2.open ("positions.m");
    for(int i=1 ; i<=nP ; i++)
      {
        outfile2 << "Pos( " << i<< ",:) = [ " ;
        outfile2 << hostCudaParticleVector[i-1].x << " " << hostCudaParticleVector[i-1].y << " " << hostCudaParticleVector[i-1].z << " ];" << std::endl;
      }
       outfile2.close();
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

    cpu_times copyToHostTime = timer.elapsed();

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "Copy to host, bin and output time: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';
    std::cout << "Total ODE integration time: " << moveTime*1e-9 << '\n';
    std::cout << "Total geometry checking time: " << geomCheckTime*1e-9 << '\n';
    std::cout << "Total ionization time: " << ionizTime*1e-9 << '\n';
    return 0;
}
