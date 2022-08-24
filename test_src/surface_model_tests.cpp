/* Sample the out-distribution? */
/* bin the in-distribution? */
#include <iostream>
#include <fstream>
#include <thrust/execution_policy.h>
#include "catch2/catch_all.hpp"
#include "config_interface.h"
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "boris.h"
#include "Surfaces.h"
#include "geometryCheck.h"
#include "boundaryInit.h"
#include "curandInitialize.h"
#include "surfaceModel.h"

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T epsilon, T margin)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) 
  {
    
    bool margin_check = (a[i] != Approx(b[i]).margin(margin));
    bool epsilon_check = (a[i] != Approx(b[i]).epsilon(epsilon));
    
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

TEST_CASE( "surface model" )
{
  SECTION( "surface model" )
  {
    /* timesteps */
    int nT = 1;
    
    /* create particles */
    libconfig::Config cfg;

    cfg.setAutoConvert(true);

    importLibConfig(cfg, SURFACE_MODEL_TEST_FILE );

    auto gitr_flags = new Flags( cfg );

    libconfig::Setting &impurity = cfg.lookup( "impurityParticleSource" );

    int nParticles = 1000;
    /* 2nd argument is deprecated - random number stream related */
    auto particleArray = new Particles( nParticles, 1, cfg, gitr_flags );

    /* equation to convert eV energy to vector velocity components */
    gitr_precision E = 200;
    gitr_precision amu = 99;
    gitr_precision vtotal = std::sqrt(2.0 * E * 1.602e-19 / amu / 1.66e-27);
    std::cout << "vtotal: " << vtotal << std::endl;

    /* set particle properties: */
    gitr_precision dt = 1.0e-2;
    gitr_precision atomic_number = 42;
    gitr_precision particle_charge = 1;
    for (int i=0; i< nParticles; i++)
    {
    particleArray->setParticleV( i, 0.5, 0.0, 1.0e-3, 0, 0, -vtotal, atomic_number, amu, particle_charge, dt );
    }
    /*
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt)
    */

    thrust::counting_iterator<std::size_t> particle_iterator_start(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(nParticles);

    int nLines = 1;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg, boundaries );

    REQUIRE( nSurfaces == 1 );

    /* start */
    int nHashes = 1;
    sim::Array<int> nR_closeGeom(nHashes, 0);
    sim::Array<int> nY_closeGeom(nHashes, 0);
    sim::Array<int> nZ_closeGeom(nHashes, 0);
    sim::Array<int> nHashPoints(nHashes, 0);
    sim::Array<int> n_closeGeomElements(nHashes, 0);
    int nEdist = 1;
    gitr_precision E0dist = 0.0;
    gitr_precision Edist = 0.0;
    int nAdist = 1;
    gitr_precision A0dist = 0.0;
    gitr_precision Adist = 0.0;
    
    getVariable(cfg, "surfaces.flux.nE", nEdist);
    getVariable(cfg, "surfaces.flux.E0", E0dist);
    getVariable(cfg, "surfaces.flux.E", Edist);

    getVariable(cfg, "surfaces.flux.nA", nAdist);
    getVariable(cfg, "surfaces.flux.A0", A0dist);
    getVariable(cfg, "surfaces.flux.A", Adist);

    auto surfaces = new Surfaces(nSurfaces, nEdist, nAdist);
    surfaces->setSurface(nEdist, E0dist, Edist, nAdist, A0dist, Adist);
    sim::Array<gitr_precision> closeGeomGridr(1),
      closeGeomGridy(1), closeGeomGridz(1);
    sim::Array<int> closeGeom(1, 0);
    /* end */
    /* 
       USE_ADAPTIVE_DT = 0 run with both options, should get same answer. Start with 0 
       turn off all hashing stuff
    */

    gitr_precision perpDiffusionCoeff = 0.0;
    cfg.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                          perpDiffusionCoeff);

    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    /* required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1 */
    /* create a unified setup script */
    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    /* uniform bfield */
    br[ 0 ] = 0;
    /* large bfield in teslas gives smaller gyromotion radius */
    by[ 0 ] = 0;
    bz[ 0 ] = 0;

    /* for the uniform efield, set efield to 1000 in z just make the cross product geometry */
    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

    int nR_PreSheathEfield = 1;
    int nY_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;

    int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

    sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);

    sim::Array<gitr_precision> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);
    PSEr[ 0 ] = 0;
    PSEz[ 0 ] = 0;
    /* y and t */
    PSEt[ 0 ] = 0;

    int nR_closeGeom_sheath = 1;
    int nY_closeGeom_sheath = 1;
    int nZ_closeGeom_sheath = 1;
    int n_closeGeomElements_sheath = 1;
    int nGeomHash_sheath = 1;
    sim::Array<gitr_precision> closeGeomGridr_sheath(nR_closeGeom_sheath),
      closeGeomGridy_sheath(nY_closeGeom_sheath),
      closeGeomGridz_sheath(nZ_closeGeom_sheath);
    sim::Array<int> closeGeom_sheath(nGeomHash_sheath);

    /* create boris operator */
    move_boris boris( particleArray, dt, boundaries.data(), nLines, nR_Bfield, nZ_Bfield,
                      bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(), by.data(),
                      nR_PreSheathEfield, 
                      nY_PreSheathEfield,
                      nZ_PreSheathEfield,
                      &preSheathEGridr.front(), &preSheathEGridy.front(),
                      &preSheathEGridz.front(), &PSEr.front(), &PSEz.front(), &PSEt.front(),
                      nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
                      n_closeGeomElements_sheath, closeGeomGridr_sheath.data(),
                      &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
                      &closeGeom_sheath.front(), gitr_flags );

    geometry_check geometry_check0(
        particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
        nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
        n_closeGeomElements.data(), &closeGeomGridr.front(),
        &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
        nEdist, E0dist, Edist, nAdist, A0dist, Adist);

#ifdef __CUDACC__
    typedef curandState rand_type;
#else
    typedef std::mt19937 rand_type;
#endif

    sim::Array<rand_type> state1(nParticles);
    
    thrust::for_each(thrust::device, particle_iterator_start, particle_iterator_end,
                   curandInitialize<>(&state1.front(), true));
  // Surface model import
  int nE_sputtRefCoeff = 1, nA_sputtRefCoeff = 1;
  int nE_sputtRefDistIn = 1, nA_sputtRefDistIn = 1;
  int nE_sputtRefDistOut = 1, nA_sputtRefDistOut = 1;
  int nE_sputtRefDistOutRef = 1, nDistE_surfaceModelRef = 1;
  int nDistE_surfaceModel = 1, nDistA_surfaceModel = 1;
  std::string surfaceModelCfg = "surfaceModel.";
  std::string surfaceModelFile;

  /* Captain! To make this test work, you must put the filestring in test_data.
     Currently, that needed file is nowhere to be found so you should disable this test for
     now. "input_path" is set to a dummy variable */
    std::string input_path = "";

    getVariable(cfg, surfaceModelCfg + "fileString", surfaceModelFile);
    nE_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
                                      surfaceModelCfg, "nEsputtRefCoeffString");
    nA_sputtRefCoeff = getDimFromFile(cfg, input_path + surfaceModelFile,
                                      surfaceModelCfg, "nAsputtRefCoeffString");
    nE_sputtRefDistIn =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistInString");
    nA_sputtRefDistIn =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nAsputtRefDistInString");
    nE_sputtRefDistOut =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistOutString");
    nE_sputtRefDistOutRef =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nEsputtRefDistOutStringRef");
    nA_sputtRefDistOut =
        getDimFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                       "nAsputtRefDistOutString");
    nDistE_surfaceModel =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOut;
    nDistE_surfaceModelRef =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nE_sputtRefDistOutRef;
    nDistA_surfaceModel =
        nE_sputtRefDistIn * nA_sputtRefDistIn * nA_sputtRefDistOut;
    std::cout << " got dimensions of surface model " << std::endl;
  
    sim::Array<gitr_precision> E_sputtRefCoeff(nE_sputtRefCoeff),
      A_sputtRefCoeff(nA_sputtRefCoeff), Elog_sputtRefCoeff(nE_sputtRefCoeff),
      energyDistGrid01(nE_sputtRefDistOut),
      energyDistGrid01Ref(nE_sputtRefDistOutRef),
      angleDistGrid01(nA_sputtRefDistOut),
      spyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
      rfyl_surfaceModel(nE_sputtRefCoeff * nA_sputtRefCoeff),
      E_sputtRefDistIn(nE_sputtRefDistIn), A_sputtRefDistIn(nA_sputtRefDistIn),
      Elog_sputtRefDistIn(nE_sputtRefDistIn),
      E_sputtRefDistOut(nE_sputtRefDistOut),
      E_sputtRefDistOutRef(nE_sputtRefDistOutRef),
      Aphi_sputtRefDistOut(nA_sputtRefDistOut),
      Atheta_sputtRefDistOut(nA_sputtRefDistOut),
      AphiDist_Y(nDistA_surfaceModel), AthetaDist_Y(nDistA_surfaceModel),
      EDist_Y(nDistE_surfaceModel), AphiDist_R(nDistA_surfaceModel),
      AthetaDist_R(nDistA_surfaceModel), EDist_R(nDistE_surfaceModelRef),
      AphiDist_CDF_Y(nDistA_surfaceModel),
      AthetaDist_CDF_Y(nDistA_surfaceModel), EDist_CDF_Y(nDistE_surfaceModel),
      AphiDist_CDF_R(nDistA_surfaceModel),
      AthetaDist_CDF_R(nDistA_surfaceModel),
      EDist_CDF_R(nDistE_surfaceModelRef),
      AphiDist_CDF_Y_regrid(nDistA_surfaceModel),
      AthetaDist_CDF_Y_regrid(nDistA_surfaceModel),
      EDist_CDF_Y_regrid(nDistE_surfaceModel),
      AphiDist_CDF_R_regrid(nDistA_surfaceModel),
      AthetaDist_CDF_R_regrid(nDistA_surfaceModel),
      EDist_CDF_R_regrid(nDistE_surfaceModelRef);
    
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefCoeff", E_sputtRefCoeff[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "A_sputtRefCoeff", A_sputtRefCoeff[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistIn", E_sputtRefDistIn[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "A_sputtRefDistIn", A_sputtRefDistIn[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistOut", E_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "E_sputtRefDistOutRef", E_sputtRefDistOutRef[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "Aphi_sputtRefDistOut", Aphi_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "Atheta_sputtRefDistOut", Atheta_sputtRefDistOut[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "sputtYldString", spyl_surfaceModel[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "reflYldString", rfyl_surfaceModel[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "EDist_Y", EDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AphiDist_Y", AphiDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AthetaDist_Y", AthetaDist_Y[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "EDist_R", EDist_R[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AphiDist_R", AphiDist_R[0]);
    getVarFromFile(cfg, input_path + surfaceModelFile, surfaceModelCfg,
                   "AthetaDist_R", AthetaDist_R[0]);
    
    for (int i = 0; i < nE_sputtRefCoeff; i++) {
      Elog_sputtRefCoeff[i] = log10(E_sputtRefCoeff[i]);
      std::cout << " EsputtRefCoeff and Elog " << E_sputtRefCoeff[i] << " "
                << Elog_sputtRefCoeff[i] << std::endl;
    }
    for (int i = 0; i < nE_sputtRefDistIn; i++) {
      Elog_sputtRefDistIn[i] = std::log10(E_sputtRefDistIn[i]);
    }
    for (int i = 0; i < nE_sputtRefDistOut; i++) {
      energyDistGrid01[i] = i * 1.0 / nE_sputtRefDistOut;
    }
    for (int i = 0; i < nE_sputtRefDistOutRef; i++) {
      energyDistGrid01Ref[i] = i * 1.0 / nE_sputtRefDistOutRef;
    }
    for (int i = 0; i < nA_sputtRefDistOut; i++) {
      angleDistGrid01[i] = i * 1.0 / nA_sputtRefDistOut;
    }

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
              EDist_Y.data(), EDist_CDF_Y.data());

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AphiDist_Y.data(), AphiDist_CDF_Y.data());

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_Y.data(), AthetaDist_CDF_Y.data());

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
              EDist_R.data(), EDist_CDF_R.data());

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AphiDist_R.data(), AphiDist_CDF_R.data());

    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_R.data(), AthetaDist_CDF_R.data());

    /* Captain! Looks like a repeated function call here */
    make2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
              AthetaDist_R.data(), AthetaDist_CDF_R.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AphiDist_CDF_Y.data(), AphiDist_CDF_Y_regrid.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AthetaDist_CDF_Y.data(), AthetaDist_CDF_Y_regrid.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOut,
                energyDistGrid01.data(), nE_sputtRefDistOut,
                E_sputtRefDistOut[nE_sputtRefDistOut - 1], EDist_CDF_Y.data(),
                EDist_CDF_Y_regrid.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Aphi_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AphiDist_CDF_R.data(), AphiDist_CDF_R_regrid.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nA_sputtRefDistOut,
                angleDistGrid01.data(), nA_sputtRefDistOut,
                Atheta_sputtRefDistOut[nA_sputtRefDistOut - 1],
                AthetaDist_CDF_R.data(), AthetaDist_CDF_R_regrid.data());

    regrid2dCDF(nE_sputtRefDistIn, nA_sputtRefDistIn, nE_sputtRefDistOutRef,
                energyDistGrid01Ref.data(), nE_sputtRefDistOutRef,
                E_sputtRefDistOutRef[nE_sputtRefDistOutRef - 1],
                EDist_CDF_R.data(), EDist_CDF_R_regrid.data());
  
    reflection reflection0(
      particleArray, dt, &state1.front(), nLines, &boundaries[0], surfaces,
      nE_sputtRefCoeff, nA_sputtRefCoeff, A_sputtRefCoeff.data(),
      Elog_sputtRefCoeff.data(), spyl_surfaceModel.data(),
      rfyl_surfaceModel.data(), nE_sputtRefDistOut, nE_sputtRefDistOutRef,
      nA_sputtRefDistOut, nE_sputtRefDistIn, nA_sputtRefDistIn,
      Elog_sputtRefDistIn.data(), A_sputtRefDistIn.data(),
      E_sputtRefDistOut.data(), E_sputtRefDistOutRef.data(),
      Aphi_sputtRefDistOut.data(), energyDistGrid01.data(),
      energyDistGrid01Ref.data(), angleDistGrid01.data(),
      EDist_CDF_Y_regrid.data(), AphiDist_CDF_Y_regrid.data(),
      EDist_CDF_R_regrid.data(), AphiDist_CDF_R_regrid.data(), nEdist, E0dist,
      Edist, nAdist, A0dist, Adist);

    /* time loop */
    for (int tt = 0; tt < nT; tt++)
    {

      thrust::for_each( thrust::device,
                        particle_iterator_start,
                        particle_iterator_end,
                        boris );

      thrust::for_each(thrust::device,
                       particle_iterator_start, 
                       particle_iterator_end, 
                       geometry_check0 );
      
      thrust::for_each(thrust::device,
                       particle_iterator_start, 
                       particle_iterator_end, 
                       reflection0 );
    }

    /* Captain! What should this distribution match? One of the input distributions? */
    std::ofstream myfile;
  myfile.open ("sputtered_v.txt");
    for (int i = 0; i < nParticles; i++)
    {
  myfile << particleArray->vx[i] << " " << 
            particleArray->vy[i] << " " << 
            particleArray->vz[i] <<"\n";
    }
  myfile.close();
  }
}
