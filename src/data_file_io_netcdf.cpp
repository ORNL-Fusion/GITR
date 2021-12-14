#include "data_file_io_netcdf.h"

int readFileDim_netcdf(const std::string& fileName,const std::string& varName)
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

int read_ar2Input_netcdf( std::string fileName, gitr_precision *Bfield[])
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

int read_profileNs_netcdf( std::string fileName, 
                           std::string nxName, std::string nzName,int &n_x,int &n_z )
{
  // Check input file exists

  std::ifstream file(fileName.c_str());
  if(!file.good())
  {
      std::cout<<"ERROR: Cannot file input file ... "<<fileName<< std::endl;
      exit(1);
  }

  netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);

  netCDF::NcDim nc_nx(nc.getDim(nxName));
  netCDF::NcDim nc_nz(nc.getDim(nzName));
  
  n_x = nc_nx.getSize(); 
  n_z = nc_nz.getSize(); 

  nc.close();

  return(0);
}
/* here and below */

int read_profileNsChar_netcdf(const char *fileName,const char *nxName,const char *nzName,int &n_x,int &n_z ) 
{

    // Check input file exists

    //ifstream file(fileName.c_str());
    std::ifstream file(fileName);
    if(!file.good()) {
        std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
        exit(1);
    }

    //NcFile nc(fileName.c_str(), NcFile::read);
    netCDF::NcFile nc(fileName, netCDF::NcFile::read);

    netCDF::NcDim nc_nx(nc.getDim(nxName));
    netCDF::NcDim nc_nz(nc.getDim(nzName));
    
    n_x = nc_nx.getSize(); 
    n_z = nc_nz.getSize(); 


    return(0);

}

int read_profiles_netcdf( std::string fileName, int &n_x, int &n_z,std::string gridxName, sim::Array<gitr_precision>& gridx,std::string gridzName,
          sim::Array<gitr_precision>& gridz,std::string dataName, sim::Array<gitr_precision>& data) {

    // Check input file exists

    std::ifstream file(fileName.c_str());
    if(!file.good()) {
        std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
        exit(1);
    }

    netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);

    netCDF::NcVar nc_gridx(nc.getVar(gridxName));
    netCDF::NcVar nc_gridz(nc.getVar(gridzName));

    nc_gridx.getVar(&gridx[0]);
    nc_gridz.getVar(&gridz[0]);
    netCDF::NcVar nc_ne(nc.getVar(dataName));
    nc_ne.getVar(&data[0]);
    nc.close();
    return(0);

}

int read_profile2d_netcdf( std::string fileName,std::string dataName, sim::Array<gitr_precision>& data) {
    std::cout << "reading 2d profile" << std::endl;
    //NcError err(2);
    //NcError::Behavior bb= (NcError::Behavior) 0;
    //NcError err(NcError::silent_nonfatal);
    // Check input file exists

    std::ifstream file(fileName.c_str());
    if(!file.good()) {
        std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
        exit(1);
    }

    netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);


    //NcVar nc_ne(nc.getVar(dataName));
    static const int NC_ERR = 2;
    netCDF::NcVar nc_ne;
    try{
        std::cout << "inside try " << std::endl;
        nc_ne = nc.getVar(dataName);
    if(nc_ne.isNull()){std::cout << "erororororo" << std::endl;}
        std::cout << "finishing try " << std::endl;
    }
    catch(int an)//NcException& e
    {
        
        std::cout << " problem in file: " << fileName << " variable: " << dataName << std::endl;
        //e.what();
        std::cout << " problem in file: " << fileName << " variable: " << dataName << std::endl;
        //return NC_ERR;
    }
    nc_ne.getVar(&data[0]);

    return(0);

}
int read_profile3d_netcdf( std::string fileName,std::string dataName, sim::Array<int>& data) {

    // Check input file exists

    std::ifstream file(fileName.c_str());
    if(!file.good()) {
        std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
        exit(1);
    }

    netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);


    netCDF::NcVar nc_ne(nc.getVar(dataName));
    nc_ne.getVar(&data[0]);

    return(0);

}

int read_profile1d_netcdf( std::string fileName,std::string gridxName, sim::Array<gitr_precision>& gridx) {

    // Check input file exists

    std::ifstream file(fileName.c_str());
    if(!file.good()) {
        std::cout<<"ERROR: Cannot file input file ... "<<fileName<<std::endl;
        exit(1);
    }

    netCDF::NcFile nc(fileName.c_str(), netCDF::NcFile::read);

    netCDF::NcVar nc_gridx(nc.getVar(gridxName));

    nc_gridx.getVar(&gridx[0]);


    return(0);

}
