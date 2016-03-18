// lab 2 program for math 578
// 1D explicit diffusion program

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <mpi.h>
#include "h1.h"
using namespace std;

void WORKER(int nWRs,int myID)
{
 int Niparms = 11;
 int iparms[Niparms];
 int Nparms = 24;
 double parms[Nparms];
 
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
double incr = 0;
int nParticlesPerProcessor;
int surfaceIndexY;
int surfaceIndexZ;

double ***VolumeBins;
double **SurfaceBins;
	double **SurfaceBinsCharge;
	double **SurfaceBinsEnergy;

	
	
	
MPI_Bcast( &iparms, Niparms, MPI_INTEGER, 0,MPI_COMM_WORLD);


MPI_Bcast( &parms, Nparms, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD);

	nP = iparms[0];
	nDensityChargeBins = iparms[1];
	nXv = parms[2];
	nYv = iparms[3];
	nZv = iparms[4];
	nY = iparms[5];
	nZ = iparms[6]; 
	nBackgroundSpecies = iparms[7];
	ionization_nDtPerApply=iparms[8];
	collision_nDtPerApply=iparms[9];
	 nT=iparms[10];
	

	
	x_start= parms[0]; 
	y_start= parms[1]; 
	z_start= parms[2]; 
	energy_eV_x_start= parms[3]; 
	energy_eV_y_start= parms[4]; 
	energy_eV_z_start= parms[5]; 
	impurity_amu= parms[6]; 
	impurity_Z= parms[7]; 
	xMinV= parms[8]; 
	xMaxV= parms[9]; 
	 yMin= parms[10];
	 yMax= parms[11];
	 zMin= parms[12];
	 zMax= parms[13];
	 surface_dz_dx= parms[14];
	 surface_zIntercept= parms[15];
	 connectionLength= parms[16];
	 nPtsPerGyroOrbit= parms[17];
	 Bx_in= parms[18];
	 By_in= parms[19];
	 Bz_in= parms[20];
	 perDiffusionCoeff_in= parms[21];
	 densitySOLDecayLength= parms[22];
	 tempSOLDecayLength= parms[23];
	 
				MPI_Barrier( MPI_COMM_WORLD );
				
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
						SurfaceBinsCharge[i][j] = 0;
						SurfaceBinsEnergy[i][j] = 0;
						SurfaceBins[i][j] = 0;
					}
				}
				
				dt = 1e-6/nPtsPerGyroOrbit;
				nParticlesPerProcessor = nP/(nWRs);
				
				//std::cout << "workers and nppp " << nWRs << nParticlesPerProcessor << std::endl;
				
					Particle Particles[nParticlesPerProcessor];
	INIT(myID,  nWRs,nP,Particles);
	std::cout<< "ppp: " << nParticlesPerProcessor << std::endl;


	unsigned long seed=(unsigned long)(time(NULL) + myID*nWRs);
	srand(seed);
	
	//std::cout<< "Worker: " << myID <<"timenull: " <<seed<< "mult: " << (1.0*myID)/(1.0*nWRs)*32.0 << std::endl;
	
	for(int p=0 ; p<nParticlesPerProcessor ; p++)
	{
	//std::cout << "p" << p<<std::endl;
		for(int tt = 0; tt< nT; tt++)
		{
			if (Particles[p].perpDistanceToSurface >= 0.0 && Particles[p].x > xMinV
			&& Particles[p].x < xMaxV && Particles[p].y > yMin && Particles[p].y < yMax
			&& Particles[p].z > zMin && Particles[p].z < zMax)
			{
			Particles[p].BorisMove(dt,  xMinV, xMaxV, yMin, yMax, zMin, zMax);
	
			Particles[p].Ionization(dt);
			}
			
			else
			{
					surfaceIndexY = int(floor((Particles[p].y - yMin)/(yMax - yMin)*(nY) + 0.0f));
		surfaceIndexZ = int(floor((Particles[p].z - zMin)/(zMax - zMin)*(nZ) + 0.0f));
		//std::cout << "Pos " << Particles[p].y << "  " << Particles[p].z << std::endl;
		//std::cout << "indices y z " << surfaceIndexY << "  " << surfaceIndexZ << std::endl;
		incr = incr + 1.0;
		SurfaceBins[surfaceIndexY][surfaceIndexZ] +=  1.0 ;
		//std::cout << "incr " << incr <<std::endl;
		SurfaceBinsCharge[surfaceIndexY][surfaceIndexZ] += Particles[p].Z ;
		SurfaceBinsEnergy[surfaceIndexY][surfaceIndexZ] += 0.5*Particles[p].amu*1.6737236e-27*(Particles[p].vx*Particles[p].vx +  Particles[p].vy*Particles[p].vy+ Particles[p].vz*Particles[p].vz)*1.60217662e-19;
			 break;
			 }
		}
		

		
	
	}
	
	//std::cout<< "Processor done " << myID << std::endl;

			
	MPI_Barrier( MPI_COMM_WORLD );
	
// 						 for(int j=0 ; j<nY ; j++)
// 				{	std::cout<< std::endl;
// 					for(int k=0 ; k<nZ ; k++)
// 					{
// 					std::cout << "  " << SurfaceBins[j][k] ;
// 					}
// 					}
// 	
	SEND_2doutput_MPI(myID,nY, nZ,SurfaceBins);
	
 		MPI_Barrier( MPI_COMM_WORLD );
 		
//  								 for(int j=0 ; j<nY ; j++)
// 				{	std::cout<< std::endl;
// 					for(int k=0 ; k<nZ ; k++)
// 					{
// 					std::cout << "  " << SurfaceBinsCharge[j][k] ;
// 					}
// 					}
		

			SEND_2doutput_MPI(myID,nY, nZ,SurfaceBinsCharge);
	
		MPI_Barrier( MPI_COMM_WORLD );
// 			SEND_2doutput_MPI(myID,nY, nZ,SurfaceBinsEnergy);
// 	
 //		MPI_Barrier( MPI_COMM_WORLD );
}
