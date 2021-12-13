#include "data_file_io.h"

int readFileDim(const std::string& fileName,const std::string& varName)
{
  netCDF::NcFile nc(fileName, netCDF::NcFile::read);

  if(nc.isNull())
  {
    std::cout << "ERROR: Failed to open " << fileName << std::endl;
  }

  netCDF::NcDim nc_nx(nc.getDim(varName));

  int n_x = nc_nx.getSize();

  nc.close();

  return n_x;
}

int read_ar2Input( std::string fileName, gitr_precision *Bfield[])
{

  // Check input file exists

  std::ifstream file(fileName.c_str());
  if(!file.good()) 
  {
      std::cout<<"ERROR: Cannot file input file ... "<<fileName<< std::endl;
      exit(1);
  }

  netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);

  netCDF::NcDim nc_nR(nc.getDim("nR"));
  netCDF::NcDim nc_nZ(nc.getDim("nZ"));
  
  int nR = nc_nR.getSize(); 
  int nZ = nc_nZ.getSize(); 

  netCDF::NcVar nc_r(nc.getVar("r"));
  netCDF::NcVar nc_z(nc.getVar("z"));

  std::vector<gitr_precision> r;
  r.resize(nR);
  nc_r.getVar(&r[0]);

  std::vector<gitr_precision> z;
  z.resize(nZ);
  nc_z.getVar(&z[0]);


  // Allocate contiguous 2D array for netcdf to work
  gitr_precision **br = new gitr_precision*[nR];
  br[0] = new gitr_precision[nR*nZ];
  for(int i=0; i<nR; i++){
      br[i] = &br[0][i*nZ];
  }


  netCDF::NcVar nc_br(nc.getVar("br"));

  nc_br.getVar(br[0]);

  for(int i=0; i<nR; i++){
      for(int j=0; j<nZ; j++){
         Bfield[i][j] = br[j][i]; 
      }
  }

  return(0);
}
