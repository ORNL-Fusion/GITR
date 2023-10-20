#include "catch2/catch_all.hpp"
#include "atomic_data_broker.h"
#include <iostream>
#include "test_data_filepath.hpp"
#include "utils.h"

#include "utils.h"
#include "Fields.h"
#include "flags.hpp"
#include "Particles.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

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

TEST_CASE("Atomic physics", "tests") {

  int const cylsymm = 0;

  SECTION("ionize - test fixed random seeds")
  {
    atomic_data_broker data_broker;

    data_broker.run();
  
    // Compare vectors to ensure reproducibility
    gitr_precision margin = 0.00001;
    gitr_precision epsilon = 0.000001;
    REQUIRE(compareVectors<gitr_precision>(data_broker.values,data_broker.values2,epsilon,margin));
  }

  SECTION("ionize - test non-fixed random seeds")
  {
    gitr_precision margin = 0.00001;
    gitr_precision epsilon = 0.000001;

    /* Captain! new code */
    atomic_data_broker data_broker;

    data_broker.run_1();
    /* end new code */

    REQUIRE(!compareVectors<gitr_precision>(data_broker.values,data_broker.values2,
                                            epsilon,margin));
  }

  SECTION("ionize - steady state")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    
    std::vector<gitr_precision> gold(20,0.0);
    gold[0] = 0.00000000000000;
    gold[1] = 0.00000000015371;
    gold[2] = 0.00000029105935;
    gold[3] = 0.00008943562195;
    gold[4] = 0.00947643526804;
    gold[5] = 0.13146921469760;
    gold[6] = 0.70242520509664;
    gold[7] = 0.14898717275460;
    gold[8] = 0.00742343038943;
    gold[9] = 0.00012814360176;

    gitr_precision margin = 500.0/nParticles;
    gitr_precision epsilon = 0.05;

    /* Captain! new code */
    atomic_data_broker data_broker;

    std::vector< gitr_precision > charge_count = data_broker.run_2();
    /* end new code */

    REQUIRE(compareVectors<gitr_precision>(charge_count,gold,epsilon,margin));
  
  }
}











