#include <iostream>
#include "catch2/catch_all.hpp"

#include "trash_gpu.h"

TEST_CASE( "" )
{
  SECTION( "" )
  {
    run_boris gpu_boris;

    gpu_boris.run();

    REQUIRE( false );
  }
}
