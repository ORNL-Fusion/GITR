#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "h1.h"
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "cudaParticle.h"
#include "boris.h"
#include "ionize.h"
#include <algorithm>
#include <boost/timer/timer.hpp>

using namespace std;
using namespace libconfig;
using namespace boost::timer;

__host__ __device__
 cudaParticle       generateParticle(double x1){
		cudaParticle p; 
               p.x = x1;
               p.y = 2.0;
               p.z = 3.0;
		
		return p;
        };

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
double perDiffusionCoeff_in;

// Background profile values used Density, temperature interpolators are 0
// or 2
double densitySOLDecayLength;
double tempSOLDecayLength;

// Background species info
int *densityChargeBins;
int *background_Z;
double *background_amu;
double *background_flow;
double *maxDensity;
double *maxTemp_eV;

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
	
	SurfaceBins = new double*[nY];
	SurfaceBinsCharge = new double*[nY];
	SurfaceBinsEnergy = new double*[nY];

 	SurfaceBins[0] = new double[nY*nZ];
 	SurfaceBinsCharge[0] = new double[nY*nZ];
 	SurfaceBinsEnergy[0] = new double[nY*nZ];
			
	for(int i=0 ; i<nY ; i++)
				{
					SurfaceBins[i] = &SurfaceBins[0][i*nZ];
					SurfaceBinsCharge[i] = &SurfaceBinsCharge[0][i*nZ];
					SurfaceBinsEnergy[i] = &SurfaceBinsEnergy[0][i*nZ];
					
				for(int j=0 ; j<nZ ; j++)
					{
						SurfaceBins[i][j] = 0;
						SurfaceBinsCharge[i][j] = 0;
						SurfaceBinsEnergy[i][j] = 0;
					}
				}
	
	//double dt;
	//double nPtsPerGyroOrbit = cfg.lookup("timeStep.nPtsPerGyroOrbit");
	//dt = 1e-6/nPtsPerGyroOrbit;

	int nP = cfg.lookup("impurityParticleSource.nP");
 	cout << "Number of particles: " << nP << endl;				
	Particle Particles[nP];
	INIT(nP,Particles, cfg);

	unsigned long seed=(unsigned long)(time(NULL));
	srand(seed);
	
	int nT = cfg.lookup("timeStep.nT");
    cout << "Number of time steps: " << nT << endl;	
    
    int surfaceIndexY;
	int surfaceIndexZ;

	thrust::host_vector<int> H(4); 
	H[0] = 14;
 	H[1] = 20;
	H[2] = 38;
	H[3] = 46;

	thrust::device_vector<int> D = H;
	D[0] = 99;
	D[1] = 88;  

	for(int i = 0; i < D.size(); i++)
	    std::cout << "D[" << i << "] = " << D[i] << std::endl;

    float dt = 1e-6;
	cudaParticle p1(x,y,z,Ex,Ey,Ez,Z,amu,dt);

    long nParticles = 1e5;
    std::cout << "nParticles: " << nParticles << std::endl;
	thrust::host_vector<cudaParticle> hostCudaParticleVector(nParticles,p1);

	//for(int i=0; i < hostCudaParticleVector.size(); i++)
	//    std::cout << hostCudaParticleVector[i].x << std::endl;

    cpu_timer timer;

	std::cout << "Initial x position GPU: " << hostCudaParticleVector[1].x << std::endl;

	thrust::device_vector<cudaParticle> deviceCudaParticleVector = hostCudaParticleVector;

    cudaThreadSynchronize();

    cpu_times copyToDeviceTime = timer.elapsed();
    std::cout << "copyToDeviceTime: " << copyToDeviceTime.wall*1e-9 << '\n';

    thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), move_boris() );

    cudaThreadSynchronize();

    cpu_times moveTimeGPU = timer.elapsed();
    std::cout << "moveTimeGPU: " << (moveTimeGPU.wall-copyToDeviceTime.wall)*1e-9 << '\n';

    //thrust::for_each(deviceCudaParticleVector.begin(), deviceCudaParticleVector.end(), ionize() );

    //cpu_times ionizeTimeGPU = timer.elapsed();
    //std::cout << "ionizeTimeGPU: " << ionizeTimeGPU.wall*1e-9 << '\n';

    thrust::host_vector<cudaParticle> hostCudaParticleVector2 = deviceCudaParticleVector;

    cudaThreadSynchronize();

    cpu_times copyToHostTime = timer.elapsed();
    std::cout << "copyToHostTime: " << (copyToHostTime.wall-moveTimeGPU.wall)*1e-9 << '\n';
	std::cout << "Final x position GPU: " << hostCudaParticleVector2.back().x << std::endl;

    // CPU

	std::vector<cudaParticle> particleVector(nParticles,p1);

    cpu_times createParticlesTimeCPU = timer.elapsed();
    std::cout << "createParticesTimeCPU: " << (createParticlesTimeCPU.wall-copyToHostTime.wall)*1e-9 << '\n';

    std::for_each( particleVector.begin(), particleVector.end(), move_boris() );

    cpu_times moveTimeCPU = timer.elapsed();
    std::cout << "moveTimeCPU: " << (moveTimeCPU.wall-createParticlesTimeCPU.wall)*1e-9 << '\n';

    std::cout << "GPU Speedup: " << (moveTimeCPU.wall-createParticlesTimeCPU.wall) / (moveTimeGPU.wall-copyToDeviceTime.wall) << '\n';
	std::cout << "Final x position CPU: " << particleVector.back().x << std::endl;
	return 0;
}
