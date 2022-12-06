#include "catch2/catch_all.hpp"
#include "coulombCollisions.h"
#include "curandInitialize.h"
#include "curandInitialize2.h"
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

    gitr_precision nu_friction = 0.0;
    gitr_precision nu_deflection = 0.0;
    gitr_precision nu_parallel = 0.0;
    gitr_precision nu_energy = 0.0;
    gitr_precision x = 0.0;
    gitr_precision y = 0.0;
    gitr_precision z = 0.0;
    gitr_precision vx = 3226.0;
    gitr_precision vy = 0.0;
    gitr_precision vz = 0.0;
    gitr_precision charge = 1.0;
    gitr_precision amu = 184.0;

    int nR_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<gitr_precision> flowVGridr(1, 0.0);
    sim::Array<gitr_precision> flowVGridz(1, 0.0);
    sim::Array<gitr_precision> flowVr(1, 0.0);
    sim::Array<gitr_precision> flowVz(1, 0.0);
    sim::Array<gitr_precision> flowVt(1, 0.0);

    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);

    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);

    gitr_precision background_Z = 1.0;
    gitr_precision background_amu = 2.0;

    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<gitr_precision> BfieldGridR(1,0.0);
    sim::Array<gitr_precision> BfieldGridZ(1,0.0);
    sim::Array<gitr_precision> BfieldR(1,0.0);
    sim::Array<gitr_precision> BfieldZ(1,0.0);
    sim::Array<gitr_precision> BfieldT(1,0.0);
    gitr_precision T_background  = 20.0;
    
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
                             BfieldT.data(), T_background, flowv_interp, cylsymm,
                             field_aligned_values );

      std::cout << "nu_friction " << nu_friction << std::endl;
      std::cout << "nu_deflection " << nu_deflection << std::endl;
      std::cout << "nu_parallel " << nu_parallel << std::endl;
      std::cout << "nu_energy " << nu_energy << std::endl;
    std::vector<gitr_precision> vals(4,0.0);
    vals[0] = nu_friction;
    vals[1] = nu_deflection;
    vals[2] = nu_parallel;
    vals[3] = nu_energy;

    std::vector<gitr_precision> gold(4,0.0);
    gold[0] = 762.588;
    gold[1] = 3023.49;
    gold[2] = 1508.46;
    gold[3] = -3006.78;

    gitr_precision margin = 0.01;
    gitr_precision epsilon = 1.0;
    REQUIRE(compareVectors<gitr_precision>(vals,gold,epsilon,margin));
  }
  
  SECTION("Frequency Evolution")
  {

    gitr_precision nu_friction = 0.0;
    gitr_precision nu_deflection = 0.0;
    gitr_precision nu_parallel = 0.0;
    gitr_precision nu_energy = 0.0;

    gitr_precision x = 0.0;
    gitr_precision y = 0.0;
    gitr_precision z = 0.0;
    
    gitr_precision vx = 0.0;
    gitr_precision vy = 0.0;
    gitr_precision vz = 3226.0;
    
    gitr_precision charge = 1.0;
    gitr_precision amu = 184.0;
    
    int nR_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<gitr_precision> flowVGridr(1, 0.0);
    sim::Array<gitr_precision> flowVGridz(1, 0.0);
    sim::Array<gitr_precision> flowVr(1, 2000.0);
    sim::Array<gitr_precision> flowVz(1, 0.0);
    sim::Array<gitr_precision> flowVt(1, 0.0);
    
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);
    
    gitr_precision background_Z = 1.0;
    gitr_precision background_amu = 2.0;
    
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<gitr_precision> BfieldGridR(1,0.0);
    sim::Array<gitr_precision> BfieldGridZ(1,0.0);
    sim::Array<gitr_precision> BfieldR(1,0.0);
    sim::Array<gitr_precision> BfieldZ(1,0.0);
    sim::Array<gitr_precision> BfieldT(1,0.0);
    gitr_precision T_background  = 20.0;
    
    int nT = 1000;
    gitr_precision dt = 1.0e-5;
    
    for(int i=0;i<nT;i++)
    {
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
                             BfieldT.data(), T_background, flowv_interp, cylsymm,
                             field_aligned_values );
     
      vx = vx + (flowVr[0] - vx)*dt*nu_friction;
    }
    
    std::vector<gitr_precision> vals(2,0.0);
    vals[0] = vx;
    vals[1] = nu_friction;

    std::vector<gitr_precision> gold(2,0.0);
    gold[1] = 762.588;
    gold[0] = 2000.0;

    gitr_precision margin = 0.01;
    gitr_precision epsilon = 1.0;
    REQUIRE(compareVectors<gitr_precision>(vals,gold,epsilon,margin));

  }
 
  SECTION("Temperature")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string file = "../test/coulomb.cfg";
    importLibConfig(cfg, COULOMB_UNIT_TEST_FILE );
    std::string input_path = "../test/";
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    //int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
    
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);

#ifdef __CUDACC__
    sim::Array<curandState> state1(nParticles);
    typedef curandState rand_type;
#else
     typedef std::mt19937 rand_type;
    sim::Array<rand_type> state1(nParticles);
#endif

    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<rand_type>(&state1.front(), 0));
    
    gitr_precision amu = 184.0;
    
    int nR_flowV = 1;
    int nY_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<gitr_precision> flowVGridr(1, 0.0);
    sim::Array<gitr_precision> flowVGridy(1, 0.0);
    sim::Array<gitr_precision> flowVGridz(1, 0.0);
    sim::Array<gitr_precision> flowVr(1, 0.0);
    sim::Array<gitr_precision> flowVz(1, 0.0);
    sim::Array<gitr_precision> flowVt(1, 0.0);
    
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
    
    gitr_precision background_Z = 1.0;
    gitr_precision background_amu = 2.0;
    
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<gitr_precision> BfieldGridR(1,0.0);
    sim::Array<gitr_precision> BfieldGridZ(1,0.0);
    sim::Array<gitr_precision> BfieldR(1,0.0);
    sim::Array<gitr_precision> BfieldZ(1,0.0);
    sim::Array<gitr_precision> BfieldT(1,0.0);
    gitr_precision T_background  = 20.0;
    
    //Adjust
    gitr_precision dt = getVariable_cfg<gitr_precision> (cfg,"timeStep.dt");
    int nT = getVariable_cfg<int> (cfg,"timeStep.nT");
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    
    coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, BfieldGridR.data(),
      &BfieldGridZ.front(), &BfieldR.front(), &BfieldZ.front(), &BfieldT.front(),gitr_flags,
      flowv_interp, cylsymm, field_aligned_values );
    
    typedef std::chrono::high_resolution_clock gitr_time;
    auto gitr_start_clock = gitr_time::now();
    
    for(int i=0; i<nT; i++)
    {
      if(i%(nT/100) == 0) std::cout << 100.0f*i/nT << " % done" << std::endl;
      thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,coulombCollisions0);
    }
    
    auto finish_clock0nc = gitr_time::now();
    typedef std::chrono::duration<gitr_precision> fsec0nc;
    fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
    printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
    std::cout << "finished ionization " << std::endl; 
  
    int n_bins = 40;
    gitr_precision v_max = 1.5e4;
    gitr_precision v_min = -1.5e4;
    gitr_precision dv = (v_max - v_min)/n_bins;
    gitr_precision mse = 0.0;
    std::vector<gitr_precision> velocity_counts(n_bins,0);
  
    for(int i=0; i< nParticles; i++)
    {
      int bin_number = std::floor((particleArray->vx[i] - v_min)/(v_max - v_min)*n_bins);	  
      if(bin_number >= 0 && bin_number < n_bins)
      {
        velocity_counts[bin_number] = velocity_counts[bin_number]+1;
      }
    }
  
    //Analytic constants
    gitr_precision B = 5.1831e-09*amu/T_background;
    for(int i=0; i<n_bins; i++)
    {
      gitr_precision this_v = v_min + (i+0.5)*dv;
      gitr_precision analytic = std::sqrt(B/3.141592653589793)*std::exp(-B*this_v*this_v);
      gitr_precision gitr = 1.0/nParticles*velocity_counts[i]/dv;
      mse = mse + std::pow(gitr-analytic,2)/n_bins;
    }
  
    gitr_precision tolerance = 1.0e-11;
    printf("mse and tol %e %e ", mse, tolerance);
    REQUIRE(mse <= tolerance);
  }
  
  SECTION("Drag")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string file = "../test/coulomb.cfg";
    importLibConfig(cfg, COULOMB_UNIT_TEST_FILE );
    std::string input_path = "../test/";
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
    
    for(int i=0; i<nParticles; i++) particleArray->vz[i] = 4580.0;
    
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);
    
#ifdef __CUDACC__
    sim::Array<curandState> state1(nParticles);
    typedef curandState rand_type;
#else
     typedef std::mt19937 rand_type;
     sim::Array<rand_type> state1(nParticles);
#endif
    
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

    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize2<rand_type>(&state1.front(),&seed.front(), &sequence.front(),&offset.front()));
    
    gitr_precision nu_friction = 0.0;
    gitr_precision nu_deflection = 0.0;
    gitr_precision nu_parallel = 0.0;
    gitr_precision nu_energy = 0.0;
    
    gitr_precision x = 0.0;
    gitr_precision y = 0.0;
    gitr_precision z = 0.0;
    
    gitr_precision vx = 0.0;
    gitr_precision vy = 2645.0;
    gitr_precision vz = 2645.0;//4580.0;
    
    gitr_precision charge = 1.0;
    gitr_precision amu = 184.0;
    
    int nR_flowV = 1;
    int nY_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<gitr_precision> flowVGridr(1, 0.0);
    sim::Array<gitr_precision> flowVGridy(1, 0.0);
    sim::Array<gitr_precision> flowVGridz(1, 0.0);
    sim::Array<gitr_precision> flowVr(1, 2000.0);
    sim::Array<gitr_precision> flowVz(1, 0.0);
    sim::Array<gitr_precision> flowVt(1, 0.0);
    
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
    
    gitr_precision background_Z = 1.0;
    gitr_precision background_amu = 2.0;
    
    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    sim::Array<gitr_precision> BfieldGridR(1,0.0);
    sim::Array<gitr_precision> BfieldGridZ(1,0.0);
    sim::Array<gitr_precision> BfieldR(1,0.0);
    sim::Array<gitr_precision> BfieldZ(1,0.0);
    sim::Array<gitr_precision> BfieldT(1,0.0);
    gitr_precision T_background  = 20.0;
    
    
    //Adjust
    gitr_precision dt = getVariable_cfg<gitr_precision> (cfg,"timeStep.dt");
    int nT = getVariable_cfg<int> (cfg,"timeStep.nT");
    
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    
    coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, BfieldGridR.data(),
      &BfieldGridZ.front(), &BfieldR.front(), &BfieldZ.front(), &BfieldT.front(),gitr_flags,
      flowv_interp, cylsymm, field_aligned_values );
    
    typedef std::chrono::high_resolution_clock gitr_time;
    auto gitr_start_clock = gitr_time::now();
    gitr_precision ave_vx = 0.0;
    std::ofstream myfile;
    myfile.open ("drag.txt");
    for(int i=0; i<nT; i++)
    {
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
                             BfieldT.data(), T_background, flowv_interp, cylsymm,
                             field_aligned_values );

      vx = vx - nu_friction*dt*(vx - flowVr[0]); 
      vy = vy - nu_friction*dt*(vy - flowVt[0]); 
      vz = vz - nu_friction*dt*(vz - flowVz[0]); 
          
      if(i%(nT/100) == 0) std::cout << 100.0f*i/nT << " % done" << std::endl;
      
      for(int j=0;j<nParticles;j++)
      {
        ave_vx = ave_vx + particleArray->vx[j]; 
      }
     
      ave_vx = ave_vx/nParticles;
      myfile  << ave_vx << " "<< vx << std::endl;
      thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,coulombCollisions0);
    }
  
    myfile.close();
    
    auto finish_clock0nc = gitr_time::now();
    typedef std::chrono::duration<gitr_precision> fsec0nc;
    fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
    printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
    std::cout << "finished ionization " << std::endl; 
    
    int n_bins = 40;
    gitr_precision v_max = 1.5e4;
    gitr_precision v_min = -1.5e4;
    gitr_precision dv = (v_max - v_min)/n_bins;
    gitr_precision mse = 0.0;
    std::vector<gitr_precision> velocity_counts(n_bins,0);
    
    for(int i=0; i< nParticles; i++)
    {
      int bin_number = std::floor((particleArray->vx[i] - v_min)/(v_max - v_min)*n_bins);	  
      if(bin_number >= 0 && bin_number < n_bins)
      {
        velocity_counts[bin_number] = velocity_counts[bin_number]+1;
      }
    }
    
    //Analytic constants
    gitr_precision B = 5.1831e-09*amu/T_background;
    std::cout << "velocity gitr analytic" << std::endl;
    
    for(int i=0; i<n_bins; i++)
    {
      gitr_precision this_v = v_min + (i+0.5)*dv;
      gitr_precision analytic = std::sqrt(B/3.141592653589793)*std::exp(-B*this_v*this_v);
      gitr_precision gitr = 1.0/nParticles*velocity_counts[i]/dv;
      mse = mse + std::pow(gitr-analytic,2)/n_bins;
    }
    
    printf("ave_vx %e \n", ave_vx);
    REQUIRE(ave_vx == Catch::Approx(flowVr[0]).margin(100.0));
  }

}
