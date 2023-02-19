#include "catch2/catch_all.hpp"
#include <iostream>
#include <netcdf>

/* Captain! Add a few things for tests in here */
TEST_CASE( "compile netcdf" )
{
  SECTION( "t0" )
  {
    netCDF::NcFile file_handle( "captain.nc", netCDF::NcFile::replace );
    REQUIRE( false );
  }
}
