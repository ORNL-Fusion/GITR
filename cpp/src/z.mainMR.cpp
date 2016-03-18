

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>
#include "h1.h"
using namespace std;

void MASTER(int nWRs, int mster)
{

char outname[] = "Deposition.m";
char outnameCharge[] = "Charge.m";
char outnameEnergy[] = "Energy.m";

//input file for GITR0.0

// Impurity particles 

int nP;
double sourceStrength;

double x_start;
double y_start;
double z_start;

double energy_eV_x_start;
double energy_eV_y_start;
double energy_eV_z_start;

double impurity_amu;
double impurity_Z;

int nDensityChargeBins;

// Volume definition

double xMinV;
double xMaxV;

// Surface definition

double yMin;
double yMax;

double zMin;
double zMax;

// Volume grid
int nXv;
int nYv;
int nZv;

// Surface grid

int nY;
int nZ;

// Surface parameterization z = dz/dx * x + b

double surface_dz_dx;
double surface_zIntercept;

double connectionLength;

// Background species info
int nBackgroundSpecies;

// Particle time stepping control

double nPtsPerGyroOrbit;
int ionization_nDtPerApply;
int collision_nDtPerApply;
int nT;

// Constant B field value - only used when BfieldInterpolator_number = 0
double Bx_in;
double By_in;
double Bz_in;

// Perp DiffusionCoeff - only used when Diffusion interpolator is = 0
double perDiffusionCoeff_in;

// Background profile values used Density, temperature interpolators are 0
// or 2
double densitySOLDecayLength;
double tempSOLDecayLength;

int *densityChargeBins;
int *background_Z;
double *background_amu;
double *background_flow;
double *maxDensity;
double *maxTemp_eV;

double dt;


int Nparms = 24;
int Niparms = 11;


int iparms[Niparms];
double parms[Nparms];



INPUT( nP,  sourceStrength, x_start, y_start, z_start, energy_eV_x_start, energy_eV_y_start,
		energy_eV_z_start, impurity_amu,  impurity_Z, nDensityChargeBins,	 xMinV, xMaxV, yMin, yMax, zMin, zMax,	 nXv,
	 nYv, nZv, nY, nZ, surface_dz_dx, surface_zIntercept, connectionLength,
	 nBackgroundSpecies, nPtsPerGyroOrbit,	 ionization_nDtPerApply, collision_nDtPerApply, nT, Bx_in,
	 By_in, Bz_in, perDiffusionCoeff_in, densitySOLDecayLength, tempSOLDecayLength	);
std::cout << "tempSOLDecayLength" << tempSOLDecayLength << std::endl;

densityChargeBins = new int[nDensityChargeBins];

background_Z = new int[nBackgroundSpecies];
background_amu = new double[nBackgroundSpecies];
background_flow = new double[nBackgroundSpecies];
maxDensity = new double[nBackgroundSpecies];
maxTemp_eV = new double[nBackgroundSpecies];

INPUT2(nDensityChargeBins,nBackgroundSpecies,densityChargeBins,background_Z,background_amu,background_flow,maxDensity,maxTemp_eV);
  std::cout << "maxTemp" << maxTemp_eV[1] << std::endl;
  
	
	iparms[0] = nP;
	iparms[1] = nDensityChargeBins;
	iparms[2] = nXv;
	iparms[3] = nYv;
	iparms[4] = nZv;
	iparms[5] = nY;
	iparms[6] = nZ;	
	iparms[7] = nBackgroundSpecies;
	iparms[8] = ionization_nDtPerApply;
	iparms[9] = collision_nDtPerApply;
	iparms[10] = nT;
	

	
	parms[0] = x_start;
	parms[1] = y_start;
	parms[2] = z_start;
	parms[3] = energy_eV_x_start;
	parms[4] = energy_eV_y_start;
	parms[5] = energy_eV_z_start;
	parms[6] = impurity_amu;
	parms[7] = impurity_Z;
	parms[8] = xMinV;
	parms[9] = xMaxV;
	parms[10] = yMin;
	parms[11] = yMax;
	parms[12] = zMin;
	parms[13] = zMax;
	parms[14] = surface_dz_dx;
	parms[15] = surface_zIntercept;
	parms[16] = connectionLength;
	parms[17] = nPtsPerGyroOrbit;
	parms[18] = Bx_in;
	parms[19] = By_in;
	parms[20] = Bz_in;
	parms[21] = perDiffusionCoeff_in;
	parms[22] = densitySOLDecayLength;
	parms[23] = tempSOLDecayLength;

	

MPI_Bcast( &iparms, Niparms, MPI_INTEGER, 0,MPI_COMM_WORLD);
MPI_Bcast( &parms, Nparms, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD);

	MPI_Barrier( MPI_COMM_WORLD );
	
	double **MSurfaceBinsGlobal;
	double **MSurfaceBins;
	double **MSurfaceBinsCharge;
	double **MSurfaceBinsChargeGlobal;
	double **MSurfaceBinsEnergy;
	
					  MSurfaceBins = new double*[nY];
					  MSurfaceBinsCharge = new double*[nY];
					  MSurfaceBinsChargeGlobal = new double*[nY];
					  MSurfaceBinsEnergy = new double*[nY];
					  MSurfaceBinsGlobal = new double*[nY];

 				 MSurfaceBins[0] = new double[nY*nZ];
 				 MSurfaceBinsCharge[0] = new double[nY*nZ];
 				 MSurfaceBinsChargeGlobal[0] = new double[nY*nZ];
 				 MSurfaceBinsEnergy[0] = new double[nY*nZ];
 				 MSurfaceBinsGlobal[0] = new double[nY*nZ];
			
			 for(int i=0 ; i<nY ; i++)
				{
					MSurfaceBins[i] = &MSurfaceBins[0][i*nZ];
					MSurfaceBinsCharge[i] = &MSurfaceBinsCharge[0][i*nZ];
					MSurfaceBinsChargeGlobal[i] = &MSurfaceBinsChargeGlobal[0][i*nZ];
					MSurfaceBinsEnergy[i] = &MSurfaceBinsEnergy[0][i*nZ];
					MSurfaceBinsGlobal[i] = &MSurfaceBinsGlobal[0][i*nZ];
					
					for(int j=0 ; j<nZ ; j++)
					{
						MSurfaceBins[i][j] = 0;
						MSurfaceBinsCharge[i][j] = 0;
						MSurfaceBinsChargeGlobal[i][j] = 0;
						MSurfaceBinsEnergy[i][j] = 0;
						MSurfaceBinsGlobal[i][j] = 0;
					}
				}
	
		MPI_Barrier( MPI_COMM_WORLD );
	

	RECV_2doutput_MPI( nWRs, nY, nZ,MSurfaceBins,MSurfaceBinsGlobal);
	
	
		MPI_Barrier( MPI_COMM_WORLD );
			
				RECV_2doutput_MPI( nWRs, nY, nZ,MSurfaceBinsCharge,MSurfaceBinsChargeGlobal);
				
	
			MPI_Barrier( MPI_COMM_WORLD );
// 				RECV_2doutput_MPI( nWRs, nY, nZ,MSurfaceBins,MSurfaceBinsEnergy);
// 	
// 			MPI_Barrier( MPI_COMM_WORLD );
			
			OUTPUT( outname,nY, nZ, MSurfaceBinsGlobal);
			OUTPUT( outnameCharge,nY, nZ, MSurfaceBinsChargeGlobal);
			OUTPUT( outnameEnergy,nY, nZ, MSurfaceBinsEnergy);
	
}