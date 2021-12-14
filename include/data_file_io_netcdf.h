#pragma once

#include <fstream>
#include <iostream>

#include "netcdf.h"
#include "ncFile.h"
#include "ncVar.h"
#include "ncDim.h"

#include "gitr_precision.h"

#include "array.h"

template <typename T>
int readFileVar_netcdf(const std::string& fileName,
                const std::string& section,const std::string& varName,T &x ) {
       std::string profiles_folder = "output/profiles";
        // Check input file exists
       std::ifstream file(fileName);
       if(!file.good()) {
         std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
         exit(1);
       }
       else{
         std::cout << "reading " <<" "  << fileName << std::endl;
       }
 
       netCDF::NcFile nc(fileName, netCDF::NcFile::read);

       if(nc.isNull()){
       std::cout << "ERROR: Failed to open " << fileName << std::endl; 
       }

       netCDF::NcVar xx;

       try{
           xx = nc.getVar(varName);
           if(xx.isNull()){std::cout << "ERROR: could not find variable "<<
              varName << " in " << fileName << std::endl;}
       }
       catch(netCDF::exceptions::NcException& e){}
       
       int numberOfDimensions = xx.getDimCount();
       if(numberOfDimensions >1)
       {
           int nTotal = 1;
           for(int j=0;j<numberOfDimensions;j++)
           {
               nTotal = nTotal*(xx.getDim(j)).getSize(); 
           }
           xx.getVar(&x);
           if(numberOfDimensions == 2)
           {
            //OUTPUT2d(profiles_folder,section+varName+".m",
            //        (xx.getDim(0)).getSize(), (xx.getDim(1)).getSize(),
            //         &x[0]);
           }
           else if(numberOfDimensions ==3)
           {
            //OUTPUT3d(profiles_folder,section+varName+".m",
            //        (xx.getDim(0)).getSize(), (xx.getDim(1)).getSize(), 
            //        (xx.getDim(2)).getSize(), &x[0]);
           }
           return nTotal;
       }
       else
       {
         netCDF::NcDim xdim;
       try{
           xdim = xx.getDim(0);
           if(xdim.isNull()){std::cout << "ERROR: could not get dimension of variable "<<
              varName << " in " << fileName << std::endl;}
       }
       catch(netCDF::exceptions::NcException& e){}
       int xlength;
       xlength = xdim.getSize();

       xx.getVar(&x);
       std::string fullName = section+varName;
       //OUTPUT1d(profiles_folder,fullName+".m", xlength, &x[0]);
       nc.close();
       return xlength;
       } 
}

/**/
int read_profileNsChar_netcdf(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z );

int read_profile2d_netcdf( std::string fileName,std::string dataName, sim::Array<gitr_precision>& data);

int read_profile1d_netcdf( std::string fileName,std::string gridxName, sim::Array<gitr_precision>& gridx);

int read_profile3d_netcdf( std::string fileName,std::string dataName, sim::Array<int>& data);


int read_profiles_netcdf( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<gitr_precision>& gridx,std::string gridzName,
                            sim::Array<gitr_precision>& gridz, std::string dataName, sim::Array<gitr_precision>& data);

/**/

int readFileDim_netcdf(const std::string& fileName,const std::string& varName);

int read_ar2Input_netcdf( std::string fileName, gitr_precision *Bfield[]);

int read_profileNs_netcdf( std::string fileName, 
                           std::string nxName, std::string nzName,int &n_x,int &n_z );
