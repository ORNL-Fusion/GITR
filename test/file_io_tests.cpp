#define CATCH_CONFIG_MAIN
//#include "file_io.hpp"
#include <iostream>
#include <libconfig.h++>
#include <stdio.h>
#include "utils.h"
#include "catch.hpp"

TEST_CASE("Factorials are computed", "[factorial]") {
  SECTION("int")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, file);
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.thing";
    int this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    REQUIRE(this_thing == 1);
  }
  SECTION("float")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, file);
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.float";
    float this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    float tol = 1e-3;
    REQUIRE_THAT(this_thing,
                     Catch::Matchers::WithinAbs(0.12345, tol));
  }
  SECTION("double")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, file);
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.float";
    double this_thing=0;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    double tol = 1e-3;
    REQUIRE_THAT(this_thing,
                     Catch::Matchers::WithinAbs(0.12345, tol));
  }
  SECTION("string")
  {
    std::cout << "starting test " << std::endl;
    typedef double P;
    libconfig::Config cfg;
    std::string file = "../test/file.cfg";
    std::cout << "about to importLibConfig " << std::endl;
    importLibConfig(cfg, file);
    std::cout << "imported cfg " << std::endl;
    const std::string var1 = "stuff.filename";
    std::string this_thing;
    getVariable(cfg, var1,this_thing);
    std::cout << "this thing values " << this_thing << std::endl;
    REQUIRE(this_thing == "netcdf_file_py.nc");
  }
}
