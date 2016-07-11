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
#include <vector>
#include "io.hpp"
#include "testRoutine.h"
#include "testRoutineCuda.h"
#include "boundaryInit.h"

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

int main()
{

Config cfg,cfg_geom;

cfg.readFile("gitrInput.cfg");
cfg_geom.readFile("gitrGeometry.cfg");

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
#endif

#if DENSITY_INTERP == 0
int nR_Dens = 1;
int nZ_Dens = 1;
std::vector<double> DensGridr(nR_Dens), DensGridz(nZ_Dens);
std::vector<double> ni(nR_Dens*nZ_Dens), ne(nR_Dens*nZ_Dens);
ni[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.ti");
ne[0] = cfg.lookup("backgroundPlasmaProfiles.Temperature.te");
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
#endif

string profilename("profiles.nc");
string densnxstring("n_x");
string densnzstring("n_z");
string densstring("ne");
string densgridxname("gridx");
string densgridzname("gridz");
int n_x;
int n_z;
int a1 = read_profileNs(profilename,densnxstring,densnzstring,n_x,n_z);
std::vector<double> gridx(n_x), gridz(n_z);
std::vector<double> dens(n_x*n_z);
int a2 = read_profiles(profilename,n_x,n_z,densgridxname, gridx,densgridzname, gridz,densstring, dens);
thrust::device_vector<double> device_dens = dens;
thrust::device_vector<double> device_gridx = gridx;
thrust::device_vector<double> device_gridz = gridz;
double interp_val1 = interp2dCombined(2.25, 0.0, -1.3,n_x,n_z, &gridx.front(), &gridz.front(), &dens.front());
std::cout << "interpolated value " << interp_val1 << std::endl;
std::vector<double> doubleVector(5,1.1);
thrust::device_vector<double> dd = doubleVector;

std::cout << "starting print loop" << std::endl;
std::for_each(doubleVector.begin(), doubleVector.end(), test_routine(2.25, 0.0, -1.3,n_x,n_z,&gridx.front(), &gridz.front(), &dens.front()) );
for (int i=0; i<5; i++)
{
        std::cout << "gridx: " << gridx[i] << std::endl;
}

double* gxptr2 = thrust::raw_pointer_cast(device_gridx.data());
double* gzptr2 = thrust::raw_pointer_cast(device_gridz.data());
double* dtptr2 = thrust::raw_pointer_cast(device_dens.data());
thrust::for_each(dd.begin(), dd.end(), test_routinecuda(2.25, 0.0, -1.3,n_x,n_z, gxptr2, gzptr2, dtptr2) );
thrust::host_vector<double> doubleVector2 = dd;
#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif
for (int i=0; i<5; i++)
{
        std::cout << "device doubleVector values: " << doubleVector2[i] << std::endl;
}


int nCS = 74;
int nTemperaturesIonize = 24;
int nDensitiesIonize = 24;
std::vector<double> rateCoeff_Ionization(nCS*nTemperaturesIonize*nDensitiesIonize);
string ADASName("ADAS_Rates_W.nc");
string IonizCoeffString("IonizationRateCoeff");
string gridTionizeName("gridTemperature_Ionization");
string gridNionizeName("gridDensity_Ionization");
std::vector<double> gridTemperature_Ionization(nTemperaturesIonize), gridDensity_Ionization(nDensitiesIonize);

int    a3 = read_profiles(ADASName, nTemperaturesIonize,nDensitiesIonize,gridTionizeName, 
        gridTemperature_Ionization,gridNionizeName,
        gridDensity_Ionization,
        IonizCoeffString,
        rateCoeff_Ionization);
        std::cout << "Coeff vector print " << 
        rateCoeff_Ionization[0*nTemperaturesIonize*nDensitiesIonize+ 1*nTemperaturesIonize+ 0] << std::endl;
double RC1 = interpRateCoeff2d ( 0, 2.25, 0.0, -1.3,n_x,n_z, &TempGridr.front(),
                      &TempGridz.front(),&te.front(),&DensGridr.front(),&DensGridz.front(), &ne.front(),nTemperaturesIonize,nDensitiesIonize,
       &gridTemperature_Ionization.front(),&gridDensity_Ionization.front(),&rateCoeff_Ionization.front() );
std::cout << "Interpolated RC " << RC1 << std::endl;


char outname[] = "Deposition.m";
char outnameCharge[] = "Charge.m";
char outnameEnergy[] = "Energy.m";

//Geometry Definition
Setting& geom = cfg_geom.lookup("geom");
int nLines = geom["x1"].getLength();
std::cout << "Number of Geometric Objects Loaded: " << nLines << std::endl;
#ifdef __CUDACC__
        thrust::host_vector<Boundary> hostBoundaryVector(nLines+1);
#else
        std::vector<Boundary> hostBoundaryVector(nLines+1);
#endif
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
    }   
hostBoundaryVector[nLines].Z = geom["Z"][nLines];
hostBoundaryVector[nLines].y1 = geom["y1"];
hostBoundaryVector[nLines].y2 = geom["y2"];
hostBoundaryVector[nLines].periodic = geom["periodic"];
std::for_each(hostBoundaryVector.begin(), hostBoundaryVector.end()-1, boundary_init(nR_Dens,nZ_Dens,&gridx.front(),&gridz.front(),&dens.front(),nR_Bfield,nZ_Bfield,&bfieldGridr.front(),&bfieldGridz.front(),&br.front(),&bz.front(), &bt.front()) );
std::cout << "exited bound_init" << std::endl;

#ifdef __CUDACC__
    thrust::device_vector<Boundary> deviceBoundaryVector = hostBoundaryVector;
    Boundary * BoundaryDevicePointer = thrust::raw_pointer_cast(deviceBoundaryVector.data());
#else
    std::vector<Boundary> * BoundaryHostPointer = &hostBoundaryVector;    
#endif
    // Volume definition

double xMinV = cfg.lookup("volumeDefinition.xMinV");
double xMaxV = cfg.lookup("volumeDefinition.xMaxV");
    // grid
int nXv = cfg.lookup("volumeDefinition.grid.nXv");
int nYv = cfg.lookup("volumeDefinition.grid.nYv");
int nZv = cfg.lookup("volumeDefinition.grid.nZv");

// Surface definition

double yMin = cfg.lookup("surfaceDefinition.yMin");
double yMax = cfg.lookup("surfaceDefinition.yMax");

double zMin = cfg.lookup("surfaceDefinition.zMin");
double zMax  = cfg.lookup("surfaceDefinition.zMax");


// Surface parameterization z = dz/dx * x + b

double surface_dz_dx  = cfg.lookup("surfaceDefinition.planeParameterization.surface_dz_dx");
double surface_zIntercept = cfg.lookup("surfaceDefinition.planeParameterization.surface_zIntercept");

// Constant B field value - only used when BfieldInterpolator_number = 0
double Bx_in = cfg.lookup("bField.Bx_in");
double By_in = cfg.lookup("bField.By_in");
double Bz_in = cfg.lookup("bField.Bz_in");
double connectionLength = cfg.lookup("bField.connectionLength");

// Particle time stepping control

int ionization_nDtPerApply  = cfg.lookup("timeStep.ionization_nDtPerApply");
int collision_nDtPerApply  = cfg.lookup("timeStep.collision_nDtPerApply");
// Perp DiffusionCoeff - only used when Diffusion interpolator is = 0
double perDiffusionCoeff_in = cfg.lookup("perpDiffusion.perDiffusionCoeff_in");

// Background species info
int *background_Z;
double *background_amu;
double *background_flow;
double *maxDensity;
double *maxTemp_eV;

#ifdef __CUDACC__
    cout<<"Using THRUST"<<endl;
#else
    cout<<"Not using THRUST"<<endl;
#endif

Setting& backgroundPlasma = cfg.lookup("backgroundPlasma");
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

    double x = cfg.lookup("impurityParticleSource.initialConditions.x_start");
    double y = cfg.lookup("impurityParticleSource.initialConditions.y_start");
    double z = cfg.lookup("impurityParticleSource.initialConditions.z_start");
    
    double Ex = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
    double Ey = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
    double Ez = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");
    
    double amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
    double Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");
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
    Particle p1(x,y,z,Ex,Ey,Ez,Z,amu);
#endif

#ifdef __CUDACC__
      thrust::host_vector<Particle> hostCudaParticleVector(nParticles,p1);
#else
        std::vector<Particle> hostCudaParticleVector(nParticles,p1);
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
       
            cpu_timer timer;

#ifdef __CUDACC__
    thrust::device_vector<Particle> deviceCudaParticleVector = hostCudaParticleVector;
#endif

    std::uniform_real_distribution<float> dist(0,1e6);
        std::random_device rd;
        std::default_random_engine generator(rd());
    
#if USEIONIZATION > 0
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
    
    cpu_times copyToDeviceTime = timer.elapsed();
    std::cout << "Initialize rand state and copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';
    for(int tt=0; tt< nT; tt++)
    {
#ifdef __CUDACC__
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), move_boris(dt,BoundaryDevicePointer, nLines) );
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
    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), crossFieldDiffusion(dt,perDiffusionCoeff_in));
#endif
#if USECOULOMBCOLLISIONS > 0
    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), coulombCollisions(dt) );
#endif
#if USETHERMALFORCE > 0
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), thermalForce(dt) );
#endif
#else
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), move_boris(dt,hostBoundaryVector,nLines) );
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), geometry_check(nLines,hostBoundaryVector) );
#if USEIONIZATION > 0
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), ionize(dt) );
#endif
#if USERECOMBINATION > 0
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), recombine(dt) );
#endif
#if USEPERPDIFFUSION > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), crossFieldDiffusion(dt,perDiffusionCoeff_in));
#endif
#if USECOULOMBCOLLISIONS > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), coulombCollisions(dt) );
#endif
#if USETHERMALFORCE > 0
        std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), thermalForce(dt) );
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


#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif

    cpu_times copyToHostTime = timer.elapsed();

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "Copy to host, bin and output time: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';
    return 0;
}
