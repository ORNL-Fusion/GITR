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
#include <algorithm>
#include <random>
#include "Particle.h"
#include <boost/timer/timer.hpp>

#ifdef __CUDACC__
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#endif

using namespace std;
using namespace libconfig;
using namespace boost::timer;

int main()
{

Config cfg;

cfg.readFile("gitrInput.cfg");

char outname[] = "Deposition.m";
char outnameCharge[] = "Charge.m";
char outnameEnergy[] = "Energy.m";

// Volume definition

double xMinV = cfg.lookup("volumeDefinition.xMinV");
double xMaxV = cfg.lookup("volumeDefinition.xMaxV");
cout << "xMaxV  " << xMaxV << endl;
    // grid
int nXv = cfg.lookup("volumeDefinition.grid.nXv");
int nYv = cfg.lookup("volumeDefinition.grid.nYv");
int nZv = cfg.lookup("volumeDefinition.grid.nZv");

// Surface definition

double yMin = cfg.lookup("surfaceDefinition.yMin");
double yMax = cfg.lookup("surfaceDefinition.yMax");

double zMin = cfg.lookup("surfaceDefinition.zMin");
double zMax  = cfg.lookup("surfaceDefinition.zMax");

// Surface grid

int nY  = cfg.lookup("surfaceDefinition.grid.nY");
int nZ  = cfg.lookup("surfaceDefinition.grid.nZ");
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
cout << "collision_nDtPerApply  " << collision_nDtPerApply << endl;
// Perp DiffusionCoeff - only used when Diffusion interpolator is = 0
double perDiffusionCoeff_in = cfg.lookup("perpDiffusion.perDiffusionCoeff_in");

// Background profile values used Density, temperature interpolators are 0 or 2
double densitySOLDecayLength;
double tempSOLDecayLength;

// Background species info
int *densityChargeBins;
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

cout << "nS  " << nS << endl;

Setting& diagnostics = cfg.lookup("diagnostics");
int nDensityChargeBins = diagnostics["densityChargeBins"].getLength();

cout << "nDensityChargeBins  " << nDensityChargeBins << endl;

densityChargeBins = new int[nDensityChargeBins];

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

cout << maxTemp_eV[i];
 }

    double x = cfg.lookup("impurityParticleSource.initialConditions.x_start");
    double y = cfg.lookup("impurityParticleSource.initialConditions.y_start");
    double z = cfg.lookup("impurityParticleSource.initialConditions.z_start");
    
    double Ex = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_x_start");
    double Ey = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_y_start");
    double Ez = cfg.lookup("impurityParticleSource.initialConditions.energy_eV_z_start");
    
    double amu = cfg.lookup("impurityParticleSource.initialConditions.impurity_amu");
    double Z = cfg.lookup("impurityParticleSource.initialConditions.impurity_Z");

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
    
    double dt;
    double nPtsPerGyroOrbit = cfg.lookup("timeStep.nPtsPerGyroOrbit");
    dt = 1e-6/nPtsPerGyroOrbit;

    int nP = cfg.lookup("impurityParticleSource.nP");
    cout << "Number of particles: " << nP << endl;              
    long nParticles = nP;
    int nT = cfg.lookup("timeStep.nT");
    cout << "Number of time steps: " << nT << endl; 
    
    int surfaceIndexY;
    int surfaceIndexZ;

    Particle p1(x,y,z,Ex,Ey,Ez,Z,amu);

    std::cout << "nParticles: " << nParticles << std::endl;
#ifdef __CUDACC__
        thrust::host_vector<Particle> hostCudaParticleVector(nParticles,p1);
#else
        std::vector<Particle> hostCudaParticleVector(nParticles,p1);
#endif
        cpu_timer timer;

    std::cout << "Initial x position GPU: " << hostCudaParticleVector[1].x << "  " << hostCudaParticleVector[0].y << "  " << hostCudaParticleVector[0].z << "  " << hostCudaParticleVector[0].vx << "  " << hostCudaParticleVector[0].vy << "  " << hostCudaParticleVector[0].vz<< "  " << hostCudaParticleVector[0].Z << std::endl;
    
#ifdef __CUDACC__
    thrust::device_vector<Particle> deviceCudaParticleVector = hostCudaParticleVector;
#endif

    //std::uniform_real_distribution<float> dist(std::numeric_limits<float>::min(),std::numeric_limits<float>::max());
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
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), move_boris(dt) );
        thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), geometry_check(1.0) );
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
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), move_boris(dt) );
    std::for_each(hostCudaParticleVector.begin(), hostCudaParticleVector.end(), geometry_check(1.0) );
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
    std::cout << "ionizeTimeGPU: " << ionizeTimeGPU.wall*1e-9 << '\n';
#ifdef __CUDACC__
    hostCudaParticleVector = deviceCudaParticleVector;
#endif

    for(int i=0; i < hostCudaParticleVector.size(); i++){
        //std::cout << " final pos" << hostCudaParticleVector[i].x << " " << hostCudaParticleVector[i].y << " " << hostCudaParticleVector[i].z << std::endl;
        if(hostCudaParticleVector[i].hitWall == 1){
        surfaceIndexY = int(floor((hostCudaParticleVector[i].y - yMin)/(yMax - yMin)*(nY) + 0.0f));
        surfaceIndexZ = int(floor((hostCudaParticleVector[i].z - zMin)/(zMax - zMin)*(nZ) + 0.0f));
        SurfaceBins[surfaceIndexY][surfaceIndexZ] +=  1.0 ;

        SurfaceBinsCharge[surfaceIndexY][surfaceIndexZ] += hostCudaParticleVector[i].Z ;
        SurfaceBinsEnergy[surfaceIndexY][surfaceIndexZ] += 0.5*hostCudaParticleVector[i].amu*1.6737236e-27*(hostCudaParticleVector[i].vx*hostCudaParticleVector[i].vx +  hostCudaParticleVector[i].vy*hostCudaParticleVector[i].vy+ hostCudaParticleVector[i].vz*hostCudaParticleVector[i].vz)/1.60217662e-19;
        }   
    }

    OUTPUT( outname,nY, nZ, SurfaceBins);
    OUTPUT( outnameCharge,nY, nZ, SurfaceBinsCharge);
    OUTPUT( outnameEnergy,nY, nZ, SurfaceBinsEnergy);
#ifdef __CUDACC__
    cudaThreadSynchronize();
#endif

    cpu_times copyToHostTime = timer.elapsed();

    std::vector<Particle> particleVector(nParticles,p1);

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "createParticesTimeCPU: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';

    std::for_each( particleVector.begin(), particleVector.end(), move_boris(dt) );

    cpu_times moveTimeCPU = timer.elapsed();
    std::cout << "moveTimeCPU: " << (moveTimeCPU.wall-createParticlesTimeCPU.wall)*1e-9 << '\n';

    return 0;
}
