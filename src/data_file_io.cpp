/* Should just contain interface functions - wrappers to be linked against the netcdf or hdf5
   versions */
#include "data_file_io.h"

#if GITR_USE_NETCDF

#include "data_file_io_netcdf.h"

#elif GITR_USE_HDF5

#include "data_file_io_hdf5.h"

#endif

int readFileDim(const std::string& fileName,const std::string& varName)
{
  #if GITR_USE_NETCDF

  return readFileDim_netcdf( fileName, varName );

  #elif GITR_USE_HDF5

  return readFileDim_hdf5( fileName, varName );

  #endif
}


int read_ar2Input( std::string fileName, gitr_precision *Bfield[])
{
  #if GITR_USE_NETCDF

  return read_ar2Input_netcdf( fileName, Bfield );

  #elif GITR_USE_HDF5

  return read_ar2Input_hdf5( fileName, Bfield );

  #endif
}

int read_profileNs( std::string fileName, 
                    std::string nxName, std::string nzName,int &n_x,int &n_z )
{
  #if GITR_USE_NETCDF

  return read_profileNs_netcdf( fileName, nxName, nzName, n_x, n_z );

  #elif GITR_USE_HDF5

  return read_profileNs_hdf5( fileName, nxName, nzName, n_x, n_z );

  #endif
}

/* here and below */
int read_profileNsChar(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z ) 
{
  #if GITR_USE_NETCDF

  return read_profileNsChar_netcdf( fileName, nxName, nzName, n_x, n_z );

  #elif GITR_USE_HDF5


  #endif
}

int read_profiles( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<gitr_precision>& gridx,std::string gridzName,
          sim::Array<gitr_precision>& gridz,std::string dataName, sim::Array<gitr_precision>& data) 
{
  #if GITR_USE_NETCDF

  return read_profiles_netcdf(  fileName,
                             n_x, 
                             n_z,
                             gridxName,
                             gridx,
                             gridzName,
                             gridz,
                             dataName,
                             data);

  #elif GITR_USE_HDF5


  #endif
}

int read_profile2d( std::string fileName,std::string dataName, sim::Array<gitr_precision>& data) 
{
  #if GITR_USE_NETCDF

  return read_profile2d_netcdf(  fileName,
                                 dataName, 
                                 data);

  #elif GITR_USE_HDF5


  #endif
}

int read_profile3d( std::string fileName,std::string dataName, sim::Array<int>& data)
{
  #if GITR_USE_NETCDF

  return read_profile3d_netcdf( fileName, dataName,  data );

  #elif GITR_USE_HDF5


  #endif
}

int read_profile1d( std::string fileName,
                    std::string gridxName,
                    sim::Array<gitr_precision>& gridx ) 
{
  #if GITR_USE_NETCDF

  return read_profile1d_netcdf(  fileName, gridxName,  gridx );


  #elif GITR_USE_HDF5


  #endif
}
