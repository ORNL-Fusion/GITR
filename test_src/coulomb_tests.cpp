#include "coulomb_data_broker.h"
#include "catch2/catch_all.hpp"
//Captain! added utils below
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "Fields.h"
#include <thrust/execution_policy.h>
#include <fstream>
#include "test_data_filepath.hpp"
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

TEST_CASE("Coulomb collision", "tests") {

  int const flowv_interp = 0;
  int const cylsymm = 0;
  int const field_aligned_values = 0;

  SECTION("Frequency")
  {
    std::vector<gitr_precision> gold(4,0.0);
    gold[0] = 762.588;
    gold[1] = 3023.49;
    gold[2] = 1508.46;
    gold[3] = -3006.78;

    gitr_precision margin = 0.01;
    gitr_precision epsilon = 1.0;

    /* Captain! new code */
    coulomb_data_broker data_broker;

    std::vector< gitr_precision > vals = data_broker.run_1();
    /* end new code */

    REQUIRE(compareVectors<gitr_precision>(vals,gold,epsilon,margin));
  }
  
  SECTION("Frequency Evolution")
  {
    std::vector<gitr_precision> gold(2,0.0);
    gold[1] = 762.588;
    gold[0] = 2000.0;

    gitr_precision margin = 0.01;
    gitr_precision epsilon = 1.0;

    /* Captain! new code. Above is trash */
    coulomb_data_broker data_broker;
    std::vector< gitr_precision > vals = data_broker.run();
    /* end new code */
    REQUIRE(compareVectors<gitr_precision>(vals,gold,epsilon,margin));
  }
 
  SECTION("Temperature")
  {
  
    gitr_precision tolerance = 1.0e-11;

    /* Captain! new code */
    coulomb_data_broker data_broker;

    double mse = data_broker.run_2();
    /* end new code */
    printf("mse and tol %e %e ", mse, tolerance);

    REQUIRE(mse <= tolerance);
  }
  
  SECTION("Drag")
  {
    /* Captain! new code */
    coulomb_data_broker data_broker;

    double ave_vx = data_broker.run_3();
    /* end new code */
    printf("ave_vx %e \n", ave_vx);

    sim::Array<gitr_precision> flowVr(1, 2000.0);
    REQUIRE(ave_vx == Catch::Approx(flowVr[0]).margin(100.0));
  }
}




























