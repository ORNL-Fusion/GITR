#ifndef _GITRUTILSIO_
#define _GITRUTILSIO_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif
#include "netcdf.h"
#include "ncFile.h"
#include "ncVar.h"
#include "ncDim.h"
#include "Boundary.h"
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include "Particle.h"
#include "libconfig.h++"
#include <libconfig.h++>
#include "interp2d.hpp"
#include <fstream>
#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif
void  read_comand_line_args(const int argc,char** argv,int& ppn,std::string& inputFile);
void checkFlags(libconfig::Config &cfg);

void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, gitr_precision *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, gitr_precision *array2d);
void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, gitr_precision *array3d);

void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, int *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, int *array2d);
void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, int *array3d);

template <typename T>
T getVariable (libconfig::Config &cfg,const std::string& s, T &a);
template <typename T>
T getVariable_cfg (libconfig::Config &cfg,const std::string& s);

extern template int getVariable(libconfig::Config &cfg,const std::string& s, int &a);
extern template float getVariable(libconfig::Config &cfg,const std::string& s, float &a);
extern template double getVariable(libconfig::Config &cfg,const std::string& s, double &a);
extern template std::string getVariable(libconfig::Config &cfg,const std::string& s, std::string &a);

template <typename P>
P get_variable(libconfig::Config &cfg, const std::string s);

extern template int get_variable(libconfig::Config &cfg, const std::string s);
extern template float get_variable(libconfig::Config &cfg, const std::string s);
extern template double get_variable(libconfig::Config &cfg, const std::string s);
extern template const char* get_variable(libconfig::Config &cfg, const std::string s);

template <typename T>
int readFileVar(const std::string& fileName,const std::string& section,const std::string& varName,T &x ) {
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

template <typename T>
int getVarFromFile (libconfig::Config &cfg,const std::string& file,const std::string& section,
        const std::string& s, T &a)
{
  std::string str;
  getVariable(cfg,section+s,str);
  int dim = readFileVar(file,section,str,a);
  return dim;
}
int getDimFromFile(libconfig::Config &cfg,const std::string& file,const std::string& section,const std::string& s);
int make2dCDF(int nX, int nY, int nZ, gitr_precision* distribution, gitr_precision* cdf);
int regrid2dCDF(int nX, int nY, int nZ,gitr_precision* xGrid,int nNew,gitr_precision maxNew, gitr_precision* cdf, gitr_precision* cdf_regrid);
int importLibConfig(libconfig::Config &cfg,std::string filepath);
int importVectorFieldNs(libconfig::Config &cfg,std::string input_path,int interpDim,std::string fieldCfgString,int &nR, int &nY,int &nZ,std::string &fileToRead);
int importVectorField(libconfig::Config &cfg,std::string input_path,int interpDim,std::string fieldCfgString,int nR, int nY,int nZ,gitr_precision &gridR,gitr_precision &gridY,gitr_precision &gridZ,gitr_precision &r, gitr_precision &y,gitr_precision &z,std::string &fileToRead);

int importGeometry(libconfig::Config &cfg,sim::Array<Boundary> &boundaries, int use_3d_geom,
                    int cylsymm, int surface_potential );

int read_ar2Input( std::string fileName, gitr_precision *Bfield[]);

int read_profileNs( std::string fileName,std::string nzName,std::string nxName,int &n_x,int &n_z );
int read_profileNsChar(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z );
int read_profile2d( std::string fileName,std::string dataName, sim::Array<gitr_precision>& data);
int read_profile1d( std::string fileName,std::string gridxName, sim::Array<gitr_precision>& gridx);
int read_profile3d( std::string fileName,std::string dataName, sim::Array<int>& data);

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<gitr_precision>& gridx,std::string gridzName,
                            sim::Array<gitr_precision>& gridz, std::string dataName, sim::Array<gitr_precision>& data);
//void OUTPUT(char outname[],int nX, int nY, float **array2d);
//void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, float *array2d);
//void OUTPUT1d(std::string folder,std::string outname,int nX, float *array2d);
//void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, float *array3d);
//void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, int *array2d);
//void OUTPUT1d(std::string folder,std::string outname,int nX, int *array2d);
//void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, int *array3d);



int readFileDim(const std::string& fileName,const std::string& varName);
int ncdfIO(int rwCode,const std::string& fileName,std::vector< std::string> dimNames,std::vector<int> dims,
        std::vector< std::string> gridNames,std::vector<int> gridMapToDims,std::vector<gitr_precision*> pointers,
        std::vector< std::string> intVarNames,std::vector<std::vector<int>> intVarDimMap, std::vector<int*> intVarPointers);
int importHashNs(libconfig::Config &cfg,std::string input_path,int nHashes,std::string fieldCfgString,int *nR, int *nY,int *nZ,int *n,int &nRTotal,int &nYTotal,int &nZTotal,int *nHashPoints, int &nHashPointsTotal,int &nGeomHash, int use_3d_geom );
#endif
