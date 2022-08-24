#include "catch2/catch_all.hpp"
#include "ionize.h"
#include "recombine.h"
#include "utils.h"
#include <thrust/execution_policy.h>
#include "curandInitialize.h"
#include "curandInitialize2.h"
#include "Fields.h"
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

TEST_CASE("Atomic physics", "tests") {

  int const cylsymm = 0;

  SECTION("ionize - test fixed random seeds")
  {
    // Creat cfg object and set to autoconvert types
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    
    // String paths to input files - this should be replaced by imported defines
    std::string input_path = "../test_data/";
    std::string file = "ionize.cfg";
    
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = 10;
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);

#ifdef __CUDACC__
    typedef curandState rand_type;
#else
    typedef std::mt19937 rand_type;
#endif

    sim::Array<rand_type> state1(nParticles);
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<>(&state1.front(), true));
    
    // GITR background fields
    // Density = 1e19 m^{-3}
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    sim::Array<gitr_precision> ne(1, 1.0e19);
    
    // Temperature = 20 eV
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);
    
    // Ionization coefficients
    int nCS_Ionize = 1, nCS_Recombine = 1;
    std::string ionizeNDens,ionizeNTemp, recombNDens, recombNTemp;
    const char *ionizeNcs, *ionizeDensGrid,
      *ionizeTempGrid, *ionizeRCvarChar, *recombNcs,
      *recombDensGrid, *recombTempGrid, *recombRCvarChar;
    std::string ionizeFile, recombFile;
    int nTemperaturesIonize = 1, nDensitiesIonize = 1;
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
        cfg.lookupValue(
            "impurityParticleSource.recombination.nChargeStateString",
            recombNcs) &&
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
    read_profileNs( ADAS_TEST_FILE, ionizeNcs, recombNcs,
                   nCS_Ionize, nCS_Recombine);

    read_profileNs( ADAS_TEST_FILE, ionizeNDens, ionizeNTemp,
                   nDensitiesIonize, nTemperaturesIonize);

    read_profileNs( ADAS_TEST_FILE, recombNDens, recombNTemp,
                   nDensitiesRecombine, nTemperaturesRecombine);
    sim::Array<gitr_precision> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
                                         nDensitiesIonize);
    sim::Array<gitr_precision> gridTemperature_Ionization(nTemperaturesIonize),
      gridDensity_Ionization(nDensitiesIonize);
    
    read_profiles( ADAS_TEST_FILE, nTemperaturesIonize, nDensitiesIonize,
      ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
      gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
    
    // Time step for this test
    gitr_precision dt = 1.0e-5;
    //auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");

    // Managed memory vector for getting random numbers from ionization functor
    sim::Array<gitr_precision> dev_f(1,-1.0);

    // Ionize functor instance
    ionize<rand_type> ionize0(
      gitr_flags,
      particleArray,
      dt,
      &state1.front(),
      nR_Dens,
      nZ_Dens,
      &DensGridr.front(),
      &DensGridz.front(),
      &ne.front(),
      nR_Temp,
      nZ_Temp,
      &TempGridr.front(),
      &TempGridz.front(), 
      &te.front(),
      nTemperaturesIonize,
      nDensitiesIonize,
      &gridTemperature_Ionization.front(),
      &gridDensity_Ionization.front(),
      &rateCoeff_Ionization.front(),
      &dev_f.front(), cylsymm );

    // Create a vector to store values in
    std::vector<gitr_precision> values(nParticles,0.0);
    
    // Get random number for each particle and store it in values
    for (int i=0;i<nParticles;i++)
    {
        thrust::for_each(thrust::device,particle_iterator0+i, particle_iterator0+i+1,ionize0);
        gitr_precision r1 = dev_f[0];
        values[i] = r1;
        //std::cout << i << " " << r1 << std::endl;

    }
    
    // Reset random number state
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<>(&state1.front(), true));
  
    // Second vector of values to compare
    std::vector<gitr_precision> values2(nParticles,0.0);
    
    for (int i=0;i<nParticles;i++)
    {
      thrust::for_each(thrust::device,particle_iterator0+i, particle_iterator0+i+1,ionize0);
      gitr_precision r1 = dev_f[0];
      values2[i] = r1;
      //std::cout << i << " " << r1 << std::endl;
    }
  
    // Compare vectors to ensure reproducibility
    gitr_precision margin = 0.00001;
    gitr_precision epsilon = 0.000001;
    REQUIRE(compareVectors<gitr_precision>(values,values2,epsilon,margin));
  }

  SECTION("ionize - test non-fixed random seeds")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);

    std::string input_path = "../test_data/";
    std::string file = "../test_data/ionize.cfg";
    
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = 10;
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
#ifdef __CUDACC__
     typedef curandState rand_type;
#else
     typedef std::mt19937 rand_type;
#endif
    sim::Array<rand_type> state1(nParticles);
  
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);
    
    sim::Array<int> seed(nParticles,0);
    sim::Array<int> sequence(nParticles,0);
    sim::Array<int> offset(nParticles,0);
    std::random_device randDeviceInit;
    std::mt19937 s0(randDeviceInit());
    std::uniform_int_distribution<int> dist(0, 100000);
    for(int i=0;i<nParticles;i++)
    {
      int r3=dist(s0);
      int r4=dist(s0);
      int r5=dist(s0);
      sequence[i] = r3;
      offset[i] = r4;
      seed[i] = r5;
    }
    
    curandInitialize2<rand_type> rand2(&state1.front(),&seed.front(), &sequence.front(),&offset.front());

    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end, rand2);
    
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    sim::Array<gitr_precision> ne(1, 1.0e19);
    
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);
  
    int nCS_Ionize = 1, nCS_Recombine = 1;
    const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
      *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
      *recombDensGrid, *recombTempGrid, *recombRCvarChar;
    std::string ionizeFile, recombFile;
    int nTemperaturesIonize = 1, nDensitiesIonize = 1;
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
        cfg.lookupValue(
            "impurityParticleSource.recombination.nChargeStateString",
            recombNcs) &&
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
    read_profileNs(ADAS_TEST_FILE, ionizeNcs, recombNcs,
                   nCS_Ionize, nCS_Recombine);

    read_profileNs(ADAS_TEST_FILE, ionizeNDens, ionizeNTemp,
                   nDensitiesIonize, nTemperaturesIonize);

    read_profileNs(ADAS_TEST_FILE, recombNDens, recombNTemp,
                   nDensitiesRecombine, nTemperaturesRecombine);
    sim::Array<gitr_precision> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
                                         nDensitiesIonize);
    sim::Array<gitr_precision> gridTemperature_Ionization(nTemperaturesIonize),
      gridDensity_Ionization(nDensitiesIonize);
    
    read_profiles(
        ADAS_TEST_FILE, nTemperaturesIonize, nDensitiesIonize,
        ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
        gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
    
    gitr_precision dt = 1.0e-5;
    //auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");

    sim::Array<gitr_precision> dev_f(1,-1.0);
    
    ionize<rand_type> ionize0(
      gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
      &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
      &rateCoeff_Ionization.front(), &dev_f.front(), cylsymm );
  
    std::vector<gitr_precision> values(nParticles,0.0);
    for (int i=0;i<nParticles;i++)
    {
        thrust::for_each(thrust::device,particle_iterator0+i, particle_iterator0+i+1,ionize0);
        gitr_precision r1 = dev_f[0];
        values[i] = r1;
        //std::cout << i << " " << r1 << std::endl;

    }
  
    for(int i=0;i<nParticles;i++)
    {
      int r3=dist(s0);
      int r4=dist(s0);
      int r5=dist(s0);
      sequence[i] = r3;
      offset[i] = r4;
      seed[i] = r5;
    }
    

    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end, rand2);
  
    std::vector<gitr_precision> values2(nParticles,0.0);
    for (int i=0;i<nParticles;i++)
    {
        thrust::for_each(thrust::device,particle_iterator0+i, particle_iterator0+i+1,ionize0);
        gitr_precision r1 = dev_f[0];
        values2[i] = r1;
        //std::cout << i << " " << r1 << std::endl;

    }
    
    gitr_precision margin = 0.00001;
    gitr_precision epsilon = 0.000001;
    REQUIRE(!compareVectors<gitr_precision>(values,values2,epsilon,margin));
  }

  SECTION("ionize - steady state")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    
    std::string input_path = "../test_data/";
    std::string file = "ionize.cfg";
    
    importLibConfig(cfg, FIELD_UNIT_TEST_FILE_0 );
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
  
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);

#if __CUDACC__
    typedef curandState rand_type;
#else
    typedef std::mt19937 rand_type;
#endif

    sim::Array<rand_type> state1(nParticles);
    sim::Array<int> seed(nParticles,0);
    sim::Array<int> sequence(nParticles,0);
    sim::Array<int> offset(nParticles,0);
    std::random_device randDeviceInit;
    std::mt19937 s0(randDeviceInit());
    std::uniform_int_distribution<int> dist(0, 100000);
    
    for(int i=0;i<nParticles;i++)
    {
       int r3=dist(s0);
       int r4=dist(s0);
       int r5=dist(s0);
       sequence[i] = r3;
       offset[i] = r4;
       seed[i] = r5;
    }
    
    curandInitialize2<rand_type> rand2(&state1.front(),&seed.front(), &sequence.front(),&offset.front());
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,rand2);
    
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    sim::Array<gitr_precision> ne(1, 1.0e19);
    
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);

    int nCS_Ionize = 1, nCS_Recombine = 1;
    const char *ionizeNcs, *ionizeNDens, *ionizeNTemp, *ionizeDensGrid,
      *ionizeTempGrid, *ionizeRCvarChar, *recombNcs, *recombNDens, *recombNTemp,
      *recombDensGrid, *recombTempGrid, *recombRCvarChar;
    std::string ionizeFile, recombFile;
    int nTemperaturesIonize = 1, nDensitiesIonize = 1;
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
    read_profileNs(ADAS_TEST_FILE, ionizeNcs, recombNcs,
                   nCS_Ionize, nCS_Recombine);

    read_profileNs(ADAS_TEST_FILE, ionizeNDens, ionizeNTemp,
                   nDensitiesIonize, nTemperaturesIonize);

    read_profileNs(ADAS_TEST_FILE, recombNDens, recombNTemp,
                   nDensitiesRecombine, nTemperaturesRecombine);
  
    sim::Array<gitr_precision> gridTemperature_Ionization(nTemperaturesIonize),
      gridDensity_Ionization(nDensitiesIonize);
    sim::Array<gitr_precision> rateCoeff_Recombination(
      nCS_Recombine * nTemperaturesRecombine * nDensitiesRecombine);
    sim::Array<gitr_precision> gridTemperature_Recombination(nTemperaturesRecombine),
      gridDensity_Recombination(nDensitiesRecombine);
    sim::Array<gitr_precision> rateCoeff_Ionization(nCS_Ionize * nTemperaturesIonize *
                                         nDensitiesIonize);
    read_profiles(
      ADAS_TEST_FILE, nTemperaturesRecombine, nDensitiesRecombine,
      recombTempGrid, gridTemperature_Recombination, recombDensGrid,
      gridDensity_Recombination, recombRCvarChar, rateCoeff_Recombination);
    
    read_profiles(
      ADAS_TEST_FILE, nTemperaturesIonize, nDensitiesIonize,
      ionizeTempGrid, gridTemperature_Ionization, ionizeDensGrid,
      gridDensity_Ionization, ionizeRCvarChar, rateCoeff_Ionization);
    
    gitr_precision dt = getVariable_cfg<gitr_precision> (cfg,"timeStep.dt");
    int nT = getVariable_cfg<int> (cfg,"timeStep.nT");
    //auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    
    sim::Array<gitr_precision> dev_f(1,-1.0);
  
    ionize<rand_type> ionize0(
      gitr_flags,particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesIonize, nDensitiesIonize,
      &gridTemperature_Ionization.front(), &gridDensity_Ionization.front(),
      &rateCoeff_Ionization.front(),&dev_f.front(), cylsymm );
  
    recombine<rand_type> recombine0(
      particleArray, dt, &state1.front(), nR_Dens, nZ_Dens, &DensGridr.front(),
      &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp, &TempGridr.front(),
      &TempGridz.front(), &te.front(), nTemperaturesRecombine,
      nDensitiesRecombine, gridTemperature_Recombination.data(),
      gridDensity_Recombination.data(), rateCoeff_Recombination.data(),gitr_flags, cylsymm );

    typedef std::chrono::high_resolution_clock gitr_time;
    auto gitr_start_clock = gitr_time::now();
    
    for(int i=0; i<nT; i++)
    {
      if(i%(nT/10) == 0) std::cout << 100.0f*i/nT << " % done" << std::endl;
      thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,ionize0);
      thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,recombine0);
    }
    
    auto finish_clock0nc = gitr_time::now();
    typedef std::chrono::duration<gitr_precision> fsec0nc;
    fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
    printf("Time taken for ionization/recombination is %6.3f (secs) \n", fs0nc.count());
    
    std::sort(particleArray->charge.begin(), particleArray->charge.end(), std::less<int>());
    
    std::vector<gitr_precision> charge_counts(20,0);
    
    for(int i=0; i< nParticles; i++)
    {
      charge_counts[int(particleArray->charge[i])] = charge_counts[int(particleArray->charge[i])]+1.0/nParticles;
    }

    //for(int i=0; i<20; i++)
    //{
    //   std::cout << charge_counts[i] << std::endl;
    //}
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
    REQUIRE(compareVectors<gitr_precision>(charge_counts,gold,epsilon,margin));
  
  }
}
