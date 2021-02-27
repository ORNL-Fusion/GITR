#include "test_utils.hpp"
#include "coulombCollisions.h"
#include "curandInitialize.h"
#include "curandInitialize2.h"
#include "Fields.h"
#include <thrust/execution_policy.h>
#include <fstream>

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
  
  SECTION("Frequency Evolution")
  {

    float nu_friction = 0.0;
    float nu_deflection = 0.0;
    float nu_parallel = 0.0;
    float nu_energy = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 3226.0;
    float charge = 1.0;
    float amu = 184.0;
    int nR_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<float> flowVGridr(1, 0.0);
    sim::Array<float> flowVGridz(1, 0.0);
    sim::Array<float> flowVr(1, 2000.0);
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
    float v_drift = 0.0;
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
                             BfieldT.data(), T_background);
      //std::cout << "nu_friction " << nu_friction << std::endl;
      //std::cout << "nu_deflection " << nu_deflection << std::endl;
      //std::cout << "nu_parallel " << nu_parallel << std::endl;
      //std::cout << "nu_energy " << nu_energy << std::endl;
      //vx = vx*(1.0-nu_energy*dt/2.0);
      vx = vx + (flowVr[0] - vx)*dt*nu_friction;
     // std::cout << "time and new velocity " << i*dt << " " << vx << std::endl;
     std::cout << i << " " << vx << " " << (flowVr[0] - vx) << " " << nu_friction << std::endl;
    }

  }
 
  SECTION("Temperature")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string file = "../test/coulomb.cfg";
    importLibConfig(cfg, file);
    std::string input_path = "../test/";
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
    printf("vz %f ", particleArray->vz[0]);
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);
    #ifdef __CUDACC__
     //typedef curandState rand_type;
    sim::Array<curandState> state1(nParticles);
#else
     typedef std::mt19937 rand_type;
    sim::Array<rand_type> state1(nParticles);
#endif
    int example=1;
    thrust::for_each(thrust::device, particle_iterator0, particle_iterator_end,
                   curandInitialize<>(&state1.front(), 0));
    
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 3226.0;
    float charge = 1.0;
    float amu = 184.0;
    int nR_flowV = 1;
    int nY_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<float> flowVGridr(1, 0.0);
    sim::Array<float> flowVGridy(1, 0.0);
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
    
    
    //Adjust
    float dt = getVariable_cfg<float> (cfg,"timeStep.dt");
    int nT = getVariable_cfg<int> (cfg,"timeStep.nT");
    std::cout << "here 1 " << std::endl; 
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    std::cout << "here 1 " << std::endl; 
    coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, BfieldGridR.data(),
      &BfieldGridZ.front(), &BfieldR.front(), &BfieldZ.front(), &BfieldT.front());
    
    std::cout << "created functors " << std::endl; 
    std::cout << "starting collisions " << std::endl; 
typedef std::chrono::high_resolution_clock gitr_time;
auto gitr_start_clock = gitr_time::now();
  for(int i=0; i<nT; i++)
  {
	  //if(i%(nT/100) == 0) std::cout << 100.0f*i/nT << " % done" << std::endl;
    thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,coulombCollisions0);
  }
  auto finish_clock0nc = gitr_time::now();
  typedef std::chrono::duration<float> fsec0nc;
  fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
  printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
    std::cout << "finished ionization " << std::endl; 
  //for(int i=0; i< nParticles; i++)
  //{
  //  std::cout << "particle charge " << i << " " << particleArray->charge[i] << std::endl;
  //}
  
  //std::sort(particleArray->vx.begin(), particleArray->vx.end(), std::less<float>());
  int n_bins = 40;
  float v_max = 1.5e4;
  float v_min = -1.5e4;
  float dv = (v_max - v_min)/n_bins;
  float mse = 0.0;
  std::vector<float> velocity_counts(n_bins,0);
  for(int i=0; i< nParticles; i++)
  {
    int bin_number = std::floor((particleArray->vx[i] - v_min)/(v_max - v_min)*n_bins);	  
    //std::cout << "particle sorted " << i << " " << particleArray->vx[i] << std::endl;
    if(bin_number >= 0 && bin_number < n_bins){
   velocity_counts[bin_number] = velocity_counts[bin_number]+1;
    }
  }
  //Analytic constants
  float B = 5.1831e-09*amu/T_background;
  std::cout << "velocity gitr analytic" << std::endl;
  for(int i=0; i<n_bins; i++)
  {
	  float this_v = v_min + (i+0.5)*dv;
	  float analytic = std::sqrt(B/3.141592653589793)*std::exp(-B*this_v*this_v);
	  float gitr = 1.0/nParticles*velocity_counts[i]/dv;
     std::cout <<this_v << " " << gitr << " " << analytic << std::endl;
     mse = mse + std::pow(gitr-analytic,2)/n_bins;
  }
  float tolerance = 1.0e-11;
  printf("mse and tol %e %e ", mse, tolerance);
  REQUIRE(mse <= tolerance);
  }
  
  SECTION("Drag")
  {
    libconfig::Config cfg;
    cfg.setAutoConvert(true);
    std::string file = "../test/coulomb.cfg";
    importLibConfig(cfg, file);
    std::string input_path = "../test/";
  
    auto gitr_flags = new Flags(cfg);
    
    int nParticles = getVariable_cfg<int> (cfg,"impurityParticleSource.nP");
    int seed01 = getVariable_cfg<int> (cfg,"operators.ionization.seed");
    auto particleArray = new Particles(nParticles,nParticles,cfg,gitr_flags);
    for(int i=0; i<nParticles; i++) particleArray->vz[i] = 4580.0;
    printf("vx vy vz %f %f %f \n", particleArray->vx[0],particleArray->vy[0],particleArray->vz[0]);
    thrust::counting_iterator<std::size_t> particle_iterator0(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);
    #ifdef __CUDACC__
     //typedef curandState rand_type;
    sim::Array<curandState> state1(nParticles);
    #else
     typedef std::mt19937 rand_type;
    sim::Array<rand_type> state1(nParticles);
    #endif
    int example=1;
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
                   curandInitialize2<>(&state1.front(),&seed.front(), &sequence.front(),&offset.front()));
    
    float nu_friction = 0.0;
    float nu_deflection = 0.0;
    float nu_parallel = 0.0;
    float nu_energy = 0.0;
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    float vx = 0.0;
    float vy = 2645.0;
    float vz = 2645.0;//4580.0;
    float charge = 1.0;
    float amu = 184.0;
    int nR_flowV = 1;
    int nY_flowV = 1;
    int nZ_flowV = 1;
    sim::Array<float> flowVGridr(1, 0.0);
    sim::Array<float> flowVGridy(1, 0.0);
    sim::Array<float> flowVGridz(1, 0.0);
    sim::Array<float> flowVr(1, 2000.0);
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
    
    
    //Adjust
    float dt = getVariable_cfg<float> (cfg,"timeStep.dt");
    int nT = getVariable_cfg<int> (cfg,"timeStep.nT");
    std::cout << "here 1 " << std::endl; 
    auto field1 = new Field(cfg,"backgroundPlasmaProfiles.Bfield");
    std::cout << "here 1 " << std::endl; 
    coulombCollisions coulombCollisions0(
      particleArray, dt, &state1.front(), nR_flowV, nY_flowV, nZ_flowV,
      &flowVGridr.front(), &flowVGridy.front(), &flowVGridz.front(),
      &flowVr.front(), &flowVz.front(), &flowVt.front(), nR_Dens, nZ_Dens,
      &DensGridr.front(), &DensGridz.front(), &ne.front(), nR_Temp, nZ_Temp,
      &TempGridr.front(), &TempGridz.front(), ti.data(), &te.front(),
      background_Z, background_amu, nR_Bfield, nZ_Bfield, BfieldGridR.data(),
      &BfieldGridZ.front(), &BfieldR.front(), &BfieldZ.front(), &BfieldT.front());
    
    std::cout << "created functors " << std::endl; 
    std::cout << "starting collisions " << std::endl; 
    typedef std::chrono::high_resolution_clock gitr_time;
    auto gitr_start_clock = gitr_time::now();
    float ave_vx = 0.0;
    float v_total = 0.0;
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
                             BfieldT.data(), T_background);
	   vx = vx - nu_friction*dt*(vx - flowVr[0]); 
	   vy = vy - nu_friction*dt*(vy - flowVt[0]); 
	   vz = vz - nu_friction*dt*(vz - flowVz[0]); 
          //if(i%(nT/1000) == 0){
        	  //std::cout << 100.0f*i/nT << " % done" << std::endl;
        	  for(int j=0;j<nParticles;j++){
        		  ave_vx = ave_vx + particleArray->vx[j]; 
        	  }
        	  ave_vx = ave_vx/nParticles;
                  myfile  << ave_vx << " "<< vx << std::endl;
          //}
      thrust::for_each(thrust::device,particle_iterator0, particle_iterator_end,coulombCollisions0);
    }
  myfile.close();
    auto finish_clock0nc = gitr_time::now();
    typedef std::chrono::duration<float> fsec0nc;
    fsec0nc fs0nc = finish_clock0nc - gitr_start_clock;
    printf("Time taken for geometry import is %6.3f (secs) \n", fs0nc.count());
    std::cout << "finished ionization " << std::endl; 
    int n_bins = 40;
    float v_max = 1.5e4;
    float v_min = -1.5e4;
    float dv = (v_max - v_min)/n_bins;
    float mse = 0.0;
    std::vector<float> velocity_counts(n_bins,0);
    for(int i=0; i< nParticles; i++)
    {
      int bin_number = std::floor((particleArray->vx[i] - v_min)/(v_max - v_min)*n_bins);	  
      //std::cout << "particle sorted " << i << " " << particleArray->vx[i] << std::endl;
      if(bin_number >= 0 && bin_number < n_bins){
     velocity_counts[bin_number] = velocity_counts[bin_number]+1;
      }
    }
    //Analytic constants
    float B = 5.1831e-09*amu/T_background;
    std::cout << "velocity gitr analytic" << std::endl;
    for(int i=0; i<n_bins; i++)
    {
            float this_v = v_min + (i+0.5)*dv;
            float analytic = std::sqrt(B/3.141592653589793)*std::exp(-B*this_v*this_v);
            float gitr = 1.0/nParticles*velocity_counts[i]/dv;
       std::cout <<this_v << " " << gitr << " " << analytic << std::endl;
       mse = mse + std::pow(gitr-analytic,2)/n_bins;
    }
    float tolerance = 1.0e-12;
    printf("ave_vx %e \n", ave_vx);
    REQUIRE(ave_vx == Approx(flowVr[0]).margin(100.0));
  }

}
