#include <iostream>
#include "test_utils.hpp"
#include "hdf5_io.h"

TEST_CASE( "low-level hdf5 interface" )
{
  SECTION( "t0" )
  {
    std::cout << "Ahoy, Captain!" << std::endl;

    write_file( "captain.h5" );

    REQUIRE( false );
  }
}
