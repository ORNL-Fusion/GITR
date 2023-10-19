#include "catch2/catch_all.hpp"
#include "cross_field_diffusion_broker.h"
#include "test_data_filepath.hpp"
#include "Particles.h"
#include "geometryCheck.h"
#include "spectroscopy.h"
#include "config_interface.h"

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T epsilon, T margin)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) 
  {
    
    bool margin_check = (a[i] != Catch::Approx(b[i]).margin(margin));
    bool epsilon_check = (a[i] != Catch::Approx(b[i]).epsilon(epsilon));
    
    if (margin_check && epsilon_check)
    {
      
      std::cout << "margin epsilon " <<
        margin_check << " " << epsilon_check << std::endl; 
      std::cout << "Element " << i << 
        " " << a[i] << " Should == " << b[i] << std::endl;
      
      return false;
    }
  }
  
  return true;
}

TEST_CASE( "cross-field diffusion operator" )
{
  /*

    Note: for this test to work, the project must be compiled with 

    GITR_USE_PERP_DIFFUSION=1
    USE_CYLSYMM=0

  */

  int const flux_ea = 1;
  int const surface_model = 1;
  int const bfield_interp = 0;
  int const use_3d_geom = 0;
  int const geom_hash = 0;
  int const spectroscopy = 2;
  int const surface_potential = 0;

  SECTION( "straight" )
  {
    cross_field_diffusion_broker data_broker;

    double rms_error = data_broker.run();

    /* Captain! End new code */

    REQUIRE(rms_error<6.0e-4);
  }

  /*

    Note: for this test to work, the project must be compiled with 

    GITR_USE_PERP_DIFFUSION=2
    USE_CYLSYMM=1

    USE3DTETGEOM=0
    USE_SURFACE_POTENTIAL=0
    PARTICLE_SOURCE_FILE=0
    USE_ADAPTIVE_DT=0

  */
  SECTION( "curved" )
  {
    cross_field_diffusion_broker data_broker;

    double rms_error = data_broker.run_1();

    /* Captain! End new code */
    REQUIRE( rms_error < 5.0e-4 );
  }
}
