#pragma once

#include <string>
#include <fstream>

#include "gitr_precision.h"
#include "array.h"

#if GITR_USE_NETCDF

#include "data_file_io_netcdf.h"

#elif GITR_USE_HDF5

#include "data_file_io_hdf5.h"

#endif

template <typename T>
int readFileVar(const std::string& fileName,
                const std::string& section,const std::string& varName,T &x ) 
{
  #if GITR_USE_NETCDF

  return readFileVar_netcdf( fileName, section, varName, x );

  #elif GITR_USE_HDF5

  return readFileVar_hdf5( fileName, section, varName, x );

  #endif
}

int readFileDim(const std::string& fileName,const std::string& varName);

int read_ar2Input( std::string fileName, gitr_precision *Bfield[]);

int read_profileNs( std::string fileName,
                    std::string nzName,
                    std::string nxName,
                    int &n_x,
                    int &n_z );

int read_profileNsChar(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z );
int read_profile2d( std::string fileName,std::string dataName, sim::Array<gitr_precision>& data);
int read_profile1d( std::string fileName,std::string gridxName, sim::Array<gitr_precision>& gridx);
int read_profile3d( std::string fileName,std::string dataName, sim::Array<int>& data);

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<gitr_precision>& gridx,std::string gridzName,
                            sim::Array<gitr_precision>& gridz, std::string dataName, sim::Array<gitr_precision>& data);
