#include "tests_general.hpp"
#include "coulombCollisions.h"
TEST_CASE("Coulomb collision", "tests") {
  SECTION("Frequency")
  {

    float nu_friction = 0.0;
    float nu_deflection = 0.0;
    float nu_parallel = 0.0;
    float nu_energy = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    float vx = 3226.0;
    float vy = 0.0;
    float vz = 0.0;
    float charge = 1.0;
    float amu = 184.0;
    int nR_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<float> flowVGridr(1, 0.0);
    sim::Array<float> flowVGridz(1, 0.0);
    sim::Array<float> flowVr(1, 0.0);
    sim::Array<float> flowVz(1, 0.0);
    sim::Array<float> flowVt(1, 0.0);
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<float> DensGridr(1, 0.0);
    sim::Array<float> DensGridz(1, 0.0);
    sim::Array<float> ni(1, 1.0e19);
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<float> TempGridr(1, 0.0);
    sim::Array<float> TempGridz(1, 0.0);
    sim::Array<float> ti(1,20.0);
    sim::Array<float> te(1,20.0);
    float background_Z = 1.0;
    float background_amu = 2.0;
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<float> BfieldGridR(1,0.0);
    sim::Array<float> BfieldGridZ(1,0.0);
    sim::Array<float> BfieldR(1,0.0);
    sim::Array<float> BfieldZ(1,0.0);
    sim::Array<float> BfieldT(1,0.0);
    float T_background  = 20.0;
    int nT = 1000;
    float dt = 1.0e-5;
    //for(int i=0;i<nT;i++)
    //{
      getSlowDownFrequencies(nu_friction, nu_deflection, nu_parallel, nu_energy,
                             x, y, z,
                             vx, vy, vz,
                             charge, amu,
                             nR_flowV, nZ_flowV, flowVGridr.data(),
                             flowVGridz.data(), flowVr.data(),
                             flowVz.data(), flowVt.data(),
                             nR_Dens, nZ_Dens, DensGridr.data(),
                             DensGridz.data(), ni.data(), nR_Temp, nZ_Temp,
                             TempGridr.data(), TempGridz.data(), ti.data(), te.data(),
                             background_Z, background_amu,
                             nR_Bfield,
                             nZ_Bfield,
                             BfieldGridR.data(),
                             BfieldGridZ.data(),
                             BfieldR.data(),
                             BfieldZ.data(),
                             BfieldT.data(), T_background);
      std::cout << "nu_friction " << nu_friction << std::endl;
      std::cout << "nu_deflection " << nu_deflection << std::endl;
      std::cout << "nu_parallel " << nu_parallel << std::endl;
      std::cout << "nu_energy " << nu_energy << std::endl;
      vx = vx*(1.0-nu_energy*dt/2.0);
     // std::cout << "time and new velocity " << i*dt << " " << vx << std::endl;
    //}

  }
}
