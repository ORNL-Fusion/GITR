#include "h1.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cstring>
using namespace std;



//IO


void INPUT(int& nP, double& sourceStrength,double& x_start,double& y_start,double& z_start,double& energy_eV_x_start,double& energy_eV_y_start,
	double&	energy_eV_z_start,double& impurity_amu, double& impurity_Z,int& nDensityChargeBins,	
	double& xMinV,double& xMaxV,double& yMin,double& yMax,double& zMin,double& zMax,	int& nXv,
	int& nYv,int& nZv,int& nY,int& nZ,double& surface_dz_dx,double& surface_zIntercept,double& connectionLength,
	int& nBackgroundSpecies,double& nPtsPerGyroOrbit,	int& ionization_nDtPerApply,int& collision_nDtPerApply,int& nT,double& Bx_in,
	double& By_in,double& Bz_in,double& perDiffusionCoeff_in,double& densitySOLDecayLength,double& tempSOLDecayLength		)
{
string string1;

	  
	  
	  ifstream gitrInputFile;
	  gitrInputFile.open ("gitrInput.txt");

 if (gitrInputFile.is_open())
  {
	gitrInputFile >> string1;
	gitrInputFile >> nP;
		std::cout << "string" << string1 << std::endl;
			std::cout << "nP" << nP << std::endl;
	
	gitrInputFile >> string1;
	gitrInputFile >> sourceStrength;
	gitrInputFile >> string1;
	gitrInputFile >> x_start;
	gitrInputFile >> string1;
	gitrInputFile >> y_start;
	gitrInputFile >> string1;
	gitrInputFile >> z_start;
	gitrInputFile >> string1;
	gitrInputFile >> energy_eV_x_start;
	gitrInputFile >> string1;
	gitrInputFile >> energy_eV_y_start;	
	gitrInputFile >> string1;
	gitrInputFile >> energy_eV_z_start;
	gitrInputFile >> string1;
	gitrInputFile >> impurity_amu;
	gitrInputFile >> string1;
	gitrInputFile >> impurity_Z;
	gitrInputFile >> string1;
	gitrInputFile >> nDensityChargeBins;
	gitrInputFile >> string1;
	gitrInputFile >> xMinV;
	gitrInputFile >> string1;
	gitrInputFile >> xMaxV;
	gitrInputFile >> string1;
	gitrInputFile >> yMin;
	gitrInputFile >> string1;
	gitrInputFile >> yMax;
	gitrInputFile >> string1;
	gitrInputFile >> zMin;
	gitrInputFile >> string1;
	gitrInputFile >> zMax;	
	gitrInputFile >> string1;
	gitrInputFile >> nXv;
	gitrInputFile >> string1;
	gitrInputFile >> nYv;
	gitrInputFile >> string1;
	gitrInputFile >> nZv;
	gitrInputFile >> string1;
	gitrInputFile >> nY;
	gitrInputFile >> string1;
	gitrInputFile >> nZ;
	gitrInputFile >> string1;
	gitrInputFile >> surface_dz_dx;
	gitrInputFile >> string1;
	gitrInputFile >> surface_zIntercept;
	gitrInputFile >> string1;
	gitrInputFile >> connectionLength;
	gitrInputFile >> string1;
	gitrInputFile >> nBackgroundSpecies;
	gitrInputFile >> string1;
	gitrInputFile >> nPtsPerGyroOrbit;	
	gitrInputFile >> string1;
	gitrInputFile >> ionization_nDtPerApply;
	gitrInputFile >> string1;
	gitrInputFile >> collision_nDtPerApply;
	gitrInputFile >> string1;
	gitrInputFile >> nT;
	gitrInputFile >> string1;
	gitrInputFile >> Bx_in;
	gitrInputFile >> string1;
	gitrInputFile >> By_in;
	gitrInputFile >> string1;
	gitrInputFile >> Bz_in;
	gitrInputFile >> string1;
	gitrInputFile >> perDiffusionCoeff_in;
	gitrInputFile >> string1;
	gitrInputFile >> densitySOLDecayLength;
	gitrInputFile >> string1;
	gitrInputFile >> tempSOLDecayLength;		


	
    gitrInputFile.close();
  }

	  else std::cout << "Unable to open file\n"; 

}

void INPUT2(int nDensityChargeBins,int nBackgroundSpecies,int densityChargeBins[],int background_Z[], double background_amu[], double background_flow[], double maxDensity[],double maxTemp_eV[])
{
string string1;
	  ifstream gitrInputFile2;
	  gitrInputFile2.open ("gitrInput2.txt");

 if (gitrInputFile2.is_open())
  {
	gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nDensityChargeBins ; i++)
	{
	gitrInputFile2 >> densityChargeBins[i];
	}

	gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nBackgroundSpecies ; i++)
	{
	gitrInputFile2 >> background_Z[i];
	}
	
		gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nBackgroundSpecies ; i++)
	{
	gitrInputFile2 >> background_amu[i];
	}
	
			gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nBackgroundSpecies ; i++)
	{
	gitrInputFile2 >> background_flow[i];
	}
	
				gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nBackgroundSpecies ; i++)
	{
	gitrInputFile2 >> maxDensity[i];
	}
	
					gitrInputFile2 >> string1;
	
	for(int i=0 ; i<nBackgroundSpecies ; i++)
	{
	gitrInputFile2 >> maxTemp_eV[i];
	}
	
	
		std::cout << "string" << string1 << std::endl;
			std::cout << "maxTemp" << maxTemp_eV[1] << std::endl;
	
	


	
    gitrInputFile2.close();
  }

	  else std::cout << "Unable to open file\n"; 

}


//OUTPUT
void OUTPUT(char outname[],int nX, int nY, double **array2d)
{
       ofstream outfile;
				//Output


			outfile.open (outname );
			
				 for(int i=1 ; i<=nX ; i++)
				{
				outfile << "Dep( " << i<< ",:) = [ " ;
					for(int j=0 ; j<nY ; j++)
					{
					outfile << array2d[i-1][j] << "  " ;
					//std::cout << r[i] << std::endl;
					}
					outfile << "  ];" << std::endl;
				}
			outfile.close();	
		
		
}
