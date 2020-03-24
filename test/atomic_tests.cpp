#include "tests_general.hpp"
#include "ionize.h"
#include "recombine.h"
#include "utils.h"
#include "Fields.h"
#include <thrust/execution_policy.h>

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != Approx(b[i]).epsilon(0.01)) {
            std::cout << a[i] << " Should == " << b[i] << std::endl;
            return false;
        }
    }
    return true;
}

TEST_CASE("Atomic physics", "tests") {
  //SECTION("ionize - test fixed random seeds")
  //{
  //  libconfig::Config cfg;
  //  cfg.setAutoConvert(true);
  //  std::string file = "../test/ionize.cfg";
  //  importLibConfig(cfg, file);
  //  std::string input_path = "../test/";
  //
  //  auto gitr_flags = new Flags(cfg);
  //  
  //  int nParticles = 10;
  //  std::cout << " autoconvert " << cfg.getAutoConvert() << std::endl;
  //  int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
  //  auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
  //  sim::Array<std::mt19937> state1(nParticles);
  //  //std::random_device randDeviceInit;
  //  //int seed0 = 250012;
  //  //for (int i =0;i < nParticles; i++) 
  //  //{
  //  //  //std::mt19937 s0(randDeviceInit());
  //  //  std::mt19937 s0(seed0+i);
  //  //  particleArray->stream_ionize[i] = s0;
  //  //}
  //  float x = 0.0;
  //  float y = 0.0;
  //  float z = 0.0;
  //  float vx = 3226.0;
  //  float vy = 0.0;
  //  float vz = 0.0;
  //  float charge = 1.0;
  //  float amu = 184.0;
  //  int nR_flowV = 1;
  //  int nZ_flowV = 1;
  //  sim::Array<float> flowVGridr(1, 0.0);
  //  sim::Array<float> flowVGridz(1, 0.0);
  //  sim::Array<float> flowVr(1, 0.0);
  //  sim::Array<float> flowVz(1, 0.0);
  //  sim::Array<float> flowVt(1, 0.0);
  //  int nR_Dens = 1;
  //  int nZ_Dens = 1;
  //  sim::Array<float> DensGridr(1, 0.0);
  //  sim::Array<float> DensGridz(1, 0.0);
  //  sim::Array<float> ni(1, 1.0e19);
  //  sim::Array<float> ne(1, 1.0e19);
  //  int nR_Temp = 1;
  //  int nZ_Temp = 1;
  //  sim::Array<float> TempGridr(1, 0.0);
  //  sim::Array<float> TempGridz(1, 0.0);
  //  sim::Array<float> ti(1,20.0);
  //  sim::Array<float> te(1,20.0);
  //  float background_Z = 1.0;
  //  float background_amu = 2.0;
  //  int nR_Bfield = 1;
  //  int nZ_Bfield = 1;
  //  sim::Array<float> BfieldGridR(1,0.0);
  //  sim::Array<float> BfieldGridZ(1,0.0);
  //  sim::Array<float> BfieldR(1,0.0);
  //  sim::Array<float> BfieldZ(1,0.0);
  //  sim::Array<float> BfieldT(1,0.0);
  //  float T_background  = 20.0;
  //int nCS_Ionize = 1, nCS_Recombine = 1;
  //const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
  //    *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
  //    *recombDensGrid, *recombTempGrid, *recombRCvarChar;
  //std::string ionizeFile, recombFile;
  //int nTemperaturesIonize = 1, nDensitiesIonize = 1;
  //int i0, i1, i2, i3, i4;
  //int nTemperaturesRecombine = 1, nDensitiesRecombine = 1;
  //  if (cfg.lookupValue("impurityParticleSource.ionization.fileString",
  //                      ionizeFile) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.nChargeStateString",
  //                      ionizeNcs) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.DensGridString",
  //                      ionizeNDens) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.TempGridString",
  //                      ionizeNTemp) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.TempGridVarName",
  //                      ionizeTempGrid) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.DensGridVarName",
  //                      ionizeDensGrid) &&
  //      cfg.lookupValue(
  //          "impurityParticleSource.recombination.nChargeStateString",
  //          recombNcs) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.CoeffVarName",
  //                      ionizeRCvarChar)) {
  //    std::cout << "Ionization rate coefficient file: " << ionizeFile
  //              << std::endl;
  //  } else {
  //    std::cout
  //        << "ERROR: Could not get ionization string info from input file "
  //        << std::endl;
  //  }
  //  std::cout << "here 0 "<< ionizeNcs << " space " << recombNcs << " space " << nCS_Ionize << nCS_Recombine << std::endl; 
  //  i0 = read_profileNs(input_path + ionizeFile, ionizeNcs, recombNcs,
  //                      nCS_Ionize, nCS_Recombine);
  //  std::cout << "here 0.1 " << std::endl; 

  //  i1 = read_profileNs(input_path + ionizeFile, ionizeNDens, ionizeNTemp,
  //                      nDensitiesIonize, nTemperaturesIonize);
  //  std::cout << "here 0.2 " << std::endl; 

  //  //i3 = read_profileNs(input_path + recombFile, recombNDens, recombNTemp,
  //  //                    nDensitiesRecombine, nTemperaturesRecombine);
  //  std::cout << "here 0.3 " << std::endl; 
  //sim::Array<float> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
  //                                       nDensitiesIonize);
  //  std::cout << "here 0.4 " << std::endl; 
  //sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
  //    gridDensity_Ionization(nDensitiesIonize);
  //  std::cout << "here 0.5 " << std::endl; 
  //  i2 = read_profiles(
  //      input_path + ionizeFile, nTemperaturesIonize, nDensitiesIonize,
  //      ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
  //      gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
  //  std::cout << "here 0.6 " << std::endl; 
  //  float dt = 1.0e-5;
  //  std::cout << "here 1 " << std::endl; 
  //  auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
  //ionize<> ionize0(
  //    gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
  //    &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
  //    &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
  //    &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
  //    &rateCoeff_Ionization.front(),field1);
  //std::vector<float> values(nParticles,0.0);
  //for (int i=0;i<nParticles;i++)
  //{
  //      std::uniform_real_distribution<float> dist(0.0, 1.0);
  //      //FIXME: r1 needs to be gottent from the proper stream
	////float r1 = dist(particleArray->stream_ionize[i]);
	//float r1 = 0.0;
  //      values[i] = r1;
  //      std::cout << i << " " << r1 << std::endl;

  //}
  //
  //std::vector<float> values2(nParticles,0.0);
  //for (int i=0;i<nParticles;i++)
  //{
  //      std::uniform_real_distribution<float> dist(0.0, 1.0);
  //      //FIXME: r1 needs to be gottent from the proper stream
  //      //float r1 = dist(particleArray->stream_ionize[i]);
  //      float r1 = 0.0;
	//values2[i] = r1;
  //      std::cout << i << " " << r1 << std::endl;

  //}
  //REQUIRE(values == values2);
  //}
  //SECTION("ionize - test non-fixed random seeds")
  //{
  //  libconfig::Config cfg;
  //  cfg.setAutoConvert(true);
  //  std::string file = "../test/ionize.cfg";
  //  importLibConfig(cfg, file);
  //  std::string input_path = "../test/";
  //
  //  auto gitr_flags = new Flags(cfg);
  //  
  //  int nParticles = 10;
  //  std::cout << " autoconvert " << cfg.getAutoConvert() << std::endl;
  //  int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
  //  auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
  //  sim::Array<std::mt19937> state1(nParticles);
  //  //std::random_device randDeviceInit;
  //  //int seed0 = 250012;
  //  //for (int i =0;i < nParticles; i++) 
  //  //{
  //  //  //std::mt19937 s0(randDeviceInit());
  //  //  std::mt19937 s0(seed0+i);
  //  //  particleArray->stream_ionize[i] = s0;
  //  //}
  //  float x = 0.0;
  //  float y = 0.0;
  //  float z = 0.0;
  //  float vx = 3226.0;
  //  float vy = 0.0;
  //  float vz = 0.0;
  //  float charge = 1.0;
  //  float amu = 184.0;
  //  int nR_flowV = 1;
  //  int nZ_flowV = 1;
  //  sim::Array<float> flowVGridr(1, 0.0);
  //  sim::Array<float> flowVGridz(1, 0.0);
  //  sim::Array<float> flowVr(1, 0.0);
  //  sim::Array<float> flowVz(1, 0.0);
  //  sim::Array<float> flowVt(1, 0.0);
  //  int nR_Dens = 1;
  //  int nZ_Dens = 1;
  //  sim::Array<float> DensGridr(1, 0.0);
  //  sim::Array<float> DensGridz(1, 0.0);
  //  sim::Array<float> ni(1, 1.0e19);
  //  sim::Array<float> ne(1, 1.0e19);
  //  int nR_Temp = 1;
  //  int nZ_Temp = 1;
  //  sim::Array<float> TempGridr(1, 0.0);
  //  sim::Array<float> TempGridz(1, 0.0);
  //  sim::Array<float> ti(1,20.0);
  //  sim::Array<float> te(1,20.0);
  //  float background_Z = 1.0;
  //  float background_amu = 2.0;
  //  int nR_Bfield = 1;
  //  int nZ_Bfield = 1;
  //  sim::Array<float> BfieldGridR(1,0.0);
  //  sim::Array<float> BfieldGridZ(1,0.0);
  //  sim::Array<float> BfieldR(1,0.0);
  //  sim::Array<float> BfieldZ(1,0.0);
  //  sim::Array<float> BfieldT(1,0.0);
  //  float T_background  = 20.0;
  //  int nT = 1000;
  //int nCS_Ionize = 1, nCS_Recombine = 1;
  //const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
  //    *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
  //    *recombDensGrid, *recombTempGrid, *recombRCvarChar;
  //std::string ionizeFile, recombFile;
  //int nTemperaturesIonize = 1, nDensitiesIonize = 1;
  //int i0, i1, i2, i3, i4;
  //int nTemperaturesRecombine = 1, nDensitiesRecombine = 1;
  //  if (cfg.lookupValue("impurityParticleSource.ionization.fileString",
  //                      ionizeFile) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.nChargeStateString",
  //                      ionizeNcs) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.DensGridString",
  //                      ionizeNDens) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.TempGridString",
  //                      ionizeNTemp) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.TempGridVarName",
  //                      ionizeTempGrid) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.DensGridVarName",
  //                      ionizeDensGrid) &&
  //      cfg.lookupValue(
  //          "impurityParticleSource.recombination.nChargeStateString",
  //          recombNcs) &&
  //      cfg.lookupValue("impurityParticleSource.ionization.CoeffVarName",
  //                      ionizeRCvarChar)) {
  //    std::cout << "Ionization rate coefficient file: " << ionizeFile
  //              << std::endl;
  //  } else {
  //    std::cout
  //        << "ERROR: Could not get ionization string info from input file "
  //        << std::endl;
  //  }
  //  std::cout << "here 0 "<< ionizeNcs << " space " << recombNcs << " space " << nCS_Ionize << nCS_Recombine << std::endl; 
  //  i0 = read_profileNs(input_path + ionizeFile, ionizeNcs, recombNcs,
  //                      nCS_Ionize, nCS_Recombine);
  //  std::cout << "here 0.1 " << std::endl; 

  //  i1 = read_profileNs(input_path + ionizeFile, ionizeNDens, ionizeNTemp,
  //                      nDensitiesIonize, nTemperaturesIonize);
  //  std::cout << "here 0.2 " << std::endl; 

  //  //i3 = read_profileNs(input_path + recombFile, recombNDens, recombNTemp,
  //  //                    nDensitiesRecombine, nTemperaturesRecombine);
  //  std::cout << "here 0.3 " << std::endl; 
  //sim::Array<float> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
  //                                       nDensitiesIonize);
  //  std::cout << "here 0.4 " << std::endl; 
  //sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
  //    gridDensity_Ionization(nDensitiesIonize);
  //  std::cout << "here 0.5 " << std::endl; 
  //  i2 = read_profiles(
  //      input_path + ionizeFile, nTemperaturesIonize, nDensitiesIonize,
  //      ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
  //      gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
  //  std::cout << "here 0.6 " << std::endl; 
  //  float dt = 1.0e-5;
  //  std::cout << "here 1 " << std::endl; 
  //  auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
  //ionize<> ionize0(
  //    gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
  //    &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
  //    &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
  //    &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
  //    &rateCoeff_Ionization.front(),field1);
  //std::vector<float> values(nParticles,0.0);
  //for (int i=0;i<nParticles;i++)
  //{
  //      std::uniform_real_distribution<float> dist(0.0, 1.0);
  //      //FIXME: r1 needs to be gottent from the proper stream
  //      float r1 = 0.0;
  //      //float r1 = dist(particleArray->stream_ionize[i]);
  //      values[i] = r1;
  //      std::cout << i << " " << r1 << std::endl;

  //}
  //
  //std::vector<float> values2(nParticles,0.0);
  //for (int i=0;i<nParticles;i++)
  //{
  //      std::uniform_real_distribution<float> dist(0.0, 1.0);
  //      //FIXME: r1 needs to be gottent from the proper stream
  //      float r1 = 0.0;
  //      //float r1 = dist(particleArray->stream_ionize[i]);
  //      values2[i] = r1;
  //      std::cout << i << " " << r1 << std::endl;

  //}
  //REQUIRE(values != values2);
  //}
  SECTION("ionize - steady state")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string file = "../test/ionize.cfg";
    importLibConfig(cfg, file);
    std::string input_path = "../test/";
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = 100000;
    std::cout << " autoconvert " << cfg.getAutoConvert() << std::endl;
    int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
    sim::Array<std::mt19937> state1(nParticles);
    std::cout << " array mt19937 " <<  state1.size() << std::endl;
    std::random_device randDeviceInit;
    int seed0 = 250012;
    for (int i =0;i < nParticles; i++) 
    {
      std::mt19937 s0(randDeviceInit());
      //std::mt19937 s0(seed0+i);
      //particleArray->stream_ionize[i] = s0;
      state1[i] = s0;
    }
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
    sim::Array<float> ne(1, 1.0e19);
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
  int nCS_Ionize = 1, nCS_Recombine = 1;
  const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
      *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
      *recombDensGrid, *recombTempGrid, *recombRCvarChar;
  std::string ionizeFile, recombFile;
  int nTemperaturesIonize = 1, nDensitiesIonize = 1;
  int i0, i1, i2, i3, i4;
  int nTemperaturesRecombine = 1, nDensitiesRecombine = 1;
    if (cfg.lookupValue("impurityParticleSource.ionization.fileString",
                        ionizeFile) &&
        cfg.lookupValue("impurityParticleSource.ionization.nChargeStateString",
                        ionizeNcs) &&
        cfg.lookupValue("impurityParticleSource.ionization.DensGridString",
                        ionizeNDens) &&
        cfg.lookupValue("impurityParticleSource.ionization.TempGridString",
                        ionizeNTemp) &&
        cfg.lookupValue("impurityParticleSource.ionization.TempGridVarName",
                        ionizeTempGrid) &&
        cfg.lookupValue("impurityParticleSource.ionization.DensGridVarName",
                        ionizeDensGrid) &&
        cfg.lookupValue("impurityParticleSource.ionization.CoeffVarName",
                        ionizeRCvarChar)) {
      std::cout << "Ionization rate coefficient file: " << ionizeFile
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get ionization string info from input file "
          << std::endl;
    }
    if (cfg.lookupValue("impurityParticleSource.recombination.fileString",
                        recombFile) &&
        cfg.lookupValue(
            "impurityParticleSource.recombination.nChargeStateString",
            recombNcs) &&
        cfg.lookupValue("impurityParticleSource.recombination.DensGridString",
                        recombNDens) &&
        cfg.lookupValue("impurityParticleSource.recombination.TempGridString",
                        recombNTemp) &&
        cfg.lookupValue("impurityParticleSource.recombination.TempGridVarName",
                        recombTempGrid) &&
        cfg.lookupValue("impurityParticleSource.recombination.DensGridVarName",
                        recombDensGrid) &&
        cfg.lookupValue("impurityParticleSource.recombination.CoeffVarName",
                        recombRCvarChar)) {
      std::cout << "Recombination rate coefficient file: " << recombFile
                << std::endl;
    } else {
      std::cout
          << "ERROR: Could not get ionization string info from input file "
          << std::endl;
    }
    std::cout << "here 0 "<< ionizeNcs << " space " << recombNcs << " space " << nCS_Ionize << nCS_Recombine << std::endl; 
    i0 = read_profileNs(input_path + ionizeFile, ionizeNcs, recombNcs,
                        nCS_Ionize, nCS_Recombine);
    std::cout << "here 0.1 " << std::endl; 

    i1 = read_profileNs(input_path + ionizeFile, ionizeNDens, ionizeNTemp,
                        nDensitiesIonize, nTemperaturesIonize);
    std::cout << "here 0.2 " << std::endl; 

    i3 = read_profileNs(input_path + recombFile, recombNDens, recombNTemp,
                        nDensitiesRecombine, nTemperaturesRecombine);
    std::cout << "here 0.3 " << std::endl; 
  sim::Array<float> gridTemperature_Ionization(nTemperaturesIonize),
      gridDensity_Ionization(nDensitiesIonize);
  sim::Array<float> rateCoeff_Recombination(
      nCS_Recombine * nTemperaturesRecombine * nDensitiesRecombine);
  sim::Array<float> gridTemperature_Recombination(nTemperaturesRecombine),
      gridDensity_Recombination(nDensitiesRecombine);
  sim::Array<float> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
                                         nDensitiesIonize);
    i4 = read_profiles(
        input_path + recombFile, nTemperaturesRecombine, nDensitiesRecombine,
        recombTempGrid, gridTemperature_Recombination, recombDensGrid,
        gridDensity_Recombination, recombRCvarChar, rateCoeff_Recombination);
    std::cout << "here 0.4 " << std::endl; 
    std::cout << "here 0.5 " << std::endl; 
    i2 = read_profiles(
        input_path + ionizeFile, nTemperaturesIonize, nDensitiesIonize,
        ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
        gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
    std::cout << "here 0.6 " << std::endl; 
    //Adjust
    float dt = 1.0e-5;
    int nT = 1000;
    std::cout << "here 1 " << std::endl; 
    auto field1 = new Field();//cfg,"backgroundPlasmaProfiles.Bfield");
  ionize<> ionize0(
      gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
      &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
      &rateCoeff_Ionization.front(), field1);
  recombine<> recombine0(
      particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesRecombine,
      nDensitiesRecombine, gridTemperature_Recombination.data(),
      gridDensity_Recombination.data(), rateCoeff_Recombination.data());
  for(auto it=particleArray->charge.begin(); it !=particleArray->charge.end(); it++)
  {
    int it2 = int(*it);
    //std::cout << "particle charge " << *it<<" "  <<   it2 << " " << particleArray->charge[it2] << std::endl;
  }
  //std::vector<int> nothing(nParticles,0);
  std::vector<int> nothing;
  for (int i=0; i< nParticles;i++)
  {
    nothing.push_back(i);
  }
  //int count = 0;
  //for(auto it=nothing.begin(); it !=nothing.end(); it++)
  //{
  //  *it = count;
  //  count++;
  //}
  std::cout << "Starting ionize loop" << std::endl;
  typedef std::chrono::high_resolution_clock gitr_time;
  auto gitr_start_clock = gitr_time::now();
  for(int i=0; i<nT; i++)
  {
  std::cout << "i " << i << std::endl;
    thrust::for_each(thrust::host,nothing.begin(), nothing.end(),ionize0);
    thrust::for_each(thrust::host,nothing.begin(), nothing.end(),recombine0);
  }
  auto finish_clock0nc = gitr_time::now();
  typedef std::chrono::duration<float> fsec0nc;
  fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
  printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
  
  std::cout << "finished ionize loop" << std::endl;
  for(int i=0; i< nParticles; i++)
  {
    //std::cout << "particle charge " << i << " " << particleArray->charge[i] << std::endl;
  }
  
  //std::sort(particleArray->charge.begin(), particleArray->charge.end(), std::less<int>());
  std::vector<float> charge_counts(10,0);
  for(int i=0; i< nParticles; i++)
  {
    if(int(particleArray->charge[i]) < 10){ 
    //std::cout << "particle sorted " << i << " " << particleArray->charge[i] << std::endl;
   charge_counts[int(particleArray->charge[i])] = charge_counts[int(particleArray->charge[i])]+1.0/nParticles;
    }
  }

  std::cout <<  "charge 4 " << charge_counts [4] << "charge 5 " << charge_counts[5] << "charge 6 "<< charge_counts[6] << " 7 " << charge_counts[7] << " 8 " << std::endl;
  std::vector<float> gold(10,0.0);
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
  REQUIRE(Approx(charge_counts[6]) == gold[6]);
  //REQUIRE(compareVectors<float>(charge_counts,gold));
  
  }
}
