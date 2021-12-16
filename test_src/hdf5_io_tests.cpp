#include <iostream>
#include "test_utils.hpp"
#include "hdf5_io.h"

TEST_CASE( "low-level hdf5 interface" )
{
  SECTION( "t0" )
  {
    write_file( "captain.h5" );

    REQUIRE( false );
  }
}
