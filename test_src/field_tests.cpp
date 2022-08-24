#include "catch2/catch_all.hpp"
#include "ionize.h"
#include "recombine.h"
#include "utils.h"
#include "test_data_filepath.hpp"

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T scale) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if ((a[i] != Catch::Approx(b[i]).scale(0.01)) && 
	    (a[i] != Catch::Approx(b[i]).margin(2.0e-10)) &&
	    (a[i] != Catch::Approx(b[i]).epsilon(std::numeric_limits<float>::epsilon()*100))) {
            std::cout << a[i] << " Should == " << b[i] << std::endl;
            return false;
        }
    }
    return true;
}

TEST_CASE("Field struct", "tests") {
  SECTION("Field - importing")
  {
    std::string inputFile = "ionize1.cfg";
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string input_path = "../test_data/";
    //importLibConfig(cfg, input_path + inputFile);
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_1 );
    //auto gitr_flags = new Flags(cfg);
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    std::cout << " nD " << field1->nD << " " << field1->dimensions.size() << std::endl;
    std::cout << " dims " << field1->dimensions[0] << " " << field1->dimensions[1] << std::endl;
    float test1 = field1->interpolate(0.5,0.1,-0.2);
    float answer = 2*std::sqrt(0.5*0.5 + 0.1*0.1) - 0.2;
    std::cout << "test1 value " << test1 << std::endl;
    REQUIRE_THAT(test1, Catch::Matchers::WithinAbs(answer,0.00001));
  }
  
  SECTION("Field - 0D or constant interpolation")
  {
    std::string inputFile = "ionize.cfg";
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string input_path = "../test_data/";
    //importLibConfig(cfg, input_path + inputFile);
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
    //auto gitr_flags = new Flags(cfg);
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    std::cout << " nD " << field1->nD << " " << field1->dimensions.size() << std::endl;
    std::cout << " dims " << field1->dimensions[0] << " " << field1->dimensions[1] << std::endl;
    float test1 = field1->interpolate(0.5,0.1,-0.2);
    float answer = 2*std::sqrt(0.5*0.5 + 0.1*0.1) - 0.2;
    std::cout << "test1 value " << test1 << std::endl;
    REQUIRE_THAT(test1, Catch::Matchers::WithinAbs(1.234,0.00001));
  }
}
