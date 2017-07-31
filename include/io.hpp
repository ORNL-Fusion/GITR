#ifndef _IO_
#define _IO_
#include <netcdf>
#include "ncFile.h"
#include "Boundary.h"
#include <vector>
#include "libconfig.h++"
using namespace std;
using namespace netCDF;
using namespace exceptions;
using namespace netCDF::exceptions;

int importGeometry(libconfig::Config &cfg,sim::Array<Boundary> &boundaries);
int read_ar2Input( std::string fileName, float *Bfield[]);

int read_profileNs( std::string fileName,std::string nzName,std::string nxName,int &n_x,int &n_z );
int read_profileNsChar(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z );
int read_profile2d( string fileName,string dataName, sim::Array<float>& data);
int read_profile1d( string fileName,string gridxName, sim::Array<float>& gridx);
int read_profile3d( string fileName,string dataName, sim::Array<int>& data);

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<float>& gridx,std::string gridzName,
                            sim::Array<float>& gridz, std::string dataName, sim::Array<float>& data);
void OUTPUT(char outname[],int nX, int nY, float **array2d);
void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, float *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, float *array2d);
void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, float *array3d);
void OUTPUT2d(std::string folder,std::string outname,int nX, int nY, int *array2d);
void OUTPUT1d(std::string folder,std::string outname,int nX, int *array2d);
void OUTPUT3d(std::string folder,std::string outname,int nX, int nY, int nZ, int *array3d);

template <typename T>
int readFileVar(const std::string& fileName,const std::string& section,const std::string& varName,T &x ) {
       std::string profiles_folder = "profiles";
        // Check input file exists
       ifstream file(fileName);
       if(!file.good()) {
         cout<<"ERROR: Cannot file input file ... "<<fileName<<endl;
         exit(1);
       }
 
       NcFile nc(fileName, NcFile::read);

       if(nc.isNull()){
       std::cout << "ERROR: Failed to open " << fileName << std::endl; 
       }

       NcVar xx;

       try{
           xx = nc.getVar(varName);
           if(xx.isNull()){std::cout << "ERROR: could not find variable "<<
              varName << " in " << fileName << std::endl;}
       }
       catch(NcException& e){}
       
       int numberOfDimensions = xx.getDimCount();
       if(numberOfDimensions >1)
       {
           int nTotal = 1;
           for(int j=0;j<numberOfDimensions;j++)
           {
               nTotal = nTotal*(xx.getDim(j)).getSize(); 
           }
           xx.getVar(&x[0]);
           if(numberOfDimensions == 2)
           {
            OUTPUT2d(profiles_folder,section+varName+".m",
                    (xx.getDim(0)).getSize(), (xx.getDim(1)).getSize(),
                     &x.front());
           }
           else if(numberOfDimensions ==3)
           {
            OUTPUT3d(profiles_folder,section+varName+".m",
                    (xx.getDim(0)).getSize(), (xx.getDim(1)).getSize(), 
                    (xx.getDim(2)).getSize(), &x.front());
           }
           return nTotal;
       }
       else
       {
       NcDim xdim;
       try{
           xdim = xx.getDim(0);
           if(xdim.isNull()){std::cout << "ERROR: could not get dimension of variable "<<
              varName << " in " << fileName << std::endl;}
       }
       catch(NcException& e){}
       int xlength;
       xlength = xdim.getSize();

       xx.getVar(&x[0]);
       std::string fullName = section+varName;
       OUTPUT1d(profiles_folder,fullName+".m", xlength, &x.front());
       return xlength;
       } 
}


int readFileDim(const std::string& fileName,const std::string& varName);
int ncdfIO(int rwCode,const std::string& fileName,vector< std::string> dimNames,vector<int> dims,
        vector< std::string> gridNames,vector<int> gridMapToDims,vector<float*> pointers,
        vector< std::string> intVarNames,vector<vector<int>> intVarDimMap, vector<int*> intVarPointers);
#endif


