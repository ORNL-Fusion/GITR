#include "config_interface.h"

netcdf_string_query::netcdf_string_query( std::string netcdf_file_init )
  :
  netcdf_file_name( netcdf_file_init ),
  nc( netcdf_file_init, netCDF::NcFile::read )
{ 
  /* Captain! Get the netcdf stuff setup here, similar to the libconfig class */
  //NcFile nc(fileName, NcFile::read);

  if(nc.isNull())
  {
    std::cout << "Captain!!! ERROR: Failed to open " << netcdf_file_init << std::endl;
  }

  else std::cout << "Successfully opened " << netcdf_file_init << std::endl;
}

netcdf_string_query::~netcdf_string_query()
{
  nc.close();
}

libconfig_string_query::libconfig_string_query( std::string libconfig_file )
{
  /* open the file */
  try
  {
    cfg.readFile( libconfig_file.c_str() );
  }

  /* bad file */
  catch(const libconfig::FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    exit( 0 );
  }

  /* bad format */
  catch(const libconfig::ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    exit( 0 );
  }
}
