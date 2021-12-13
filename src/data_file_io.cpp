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
