#include "array.h"
#include <cmath>
#include "utils.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <netcdf.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "flags.hpp"


// Cache for surface reaction data
std::map<
        std::pair<std::string, std::string>, std::tuple<int, int, std::vector<double>, std::vector<double>, std::vector<double>
        , std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,  int, int,
        std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, 
        std::vector<double>, std::vector<double>, std::vector<double>, 
         std::vector<double>,
          std::vector<double>, std::vector<double>,        std::vector<double>, 
        std::vector<double>, std::vector<double>, int, int, int, int, int,
         std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
          std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> , gitr_precision, gitr_precision, gitr_precision, gitr_precision, gitr_precision, gitr_precision,
          int, int, gitr_precision, gitr_precision, gitr_precision, gitr_precision, int 
        //    nEdist, nAdist, E0dist, Edist, A0dist, Adist
        > > cacheSurfaceData;

inline std::tuple<int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> , std::vector<double>, std::vector<double>,  int, int,
 std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, 
 std::vector<double>, std::vector<double>, std::vector<double>,  int, int, int, int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
  std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, gitr_precision, gitr_precision, gitr_precision, gitr_precision, gitr_precision, gitr_precision,
  int, int, gitr_precision, gitr_precision, gitr_precision, gitr_precision, int 
> processSurfaceData(const std::string& incidentMaterial, const std::string& targetMaterial)
    {
        auto key = std::make_pair(incidentMaterial, targetMaterial);
        if (cacheSurfaceData.find(key) != cacheSurfaceData.end()) {
            return cacheSurfaceData[key];  // Return cached result
        }
    // Surface model import
    int nE_sputtRefCoeff = 1, nA_sputtRefCoeff = 1;
    int nE_sputtRefDistIn = 1, nA_sputtRefDistIn = 1;
    int nE_sputtRefDistOut = 1, nA_sputtRefDistOut = 1;
    int nE_sputtRefDistOutRef = 1, nDistE_surfaceModelRef = 1;
    int nDistE_surfaceModel = 1, nDistA_surfaceModel = 1;

    int nEdist = 1;
    gitr_precision E0dist = 0.0;
    gitr_precision Edist = 0.0;
    int nAdist = 1;
    gitr_precision A0dist = 0.0;
    gitr_precision Adist = 0.0;
    
    std::string surfaceModelCfg = "surfaceModel.";
    std::string surfaceModelFile;
    const std::string input_path = "input/";
    // config file is gitrInput.cfg
    libconfig::Config cfg;
    std::string cfgFile = input_path + "gitrInput.cfg";
    surfaceModelFile = "surfaceReactions_" + incidentMaterial + "_on_" + targetMaterial + ".nc";


    //cfg is the config file
    cfg.readFile(cfgFile.c_str());

    getVariable(cfg, "surfaces.flux.nE", nEdist);
    getVariable(cfg, "surfaces.flux.E0", E0dist);
    getVariable(cfg, "surfaces.flux.E", Edist);

    getVariable(cfg, "surfaces.flux.nA", nAdist);
    getVariable(cfg, "surfaces.flux.A0", A0dist);
    getVariable(cfg, "surfaces.flux.A", Adist);


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
    

    std::vector<double> E_sputtRefCoeff(nE_sputtRefCoeff),
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
    
    gitr_precision spylAInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AphiDist_CDF_Y_regrid.data());
    gitr_precision spylAthetaInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AthetaDist_CDF_Y_regrid.data());
    gitr_precision sputEInterpVal = interp3d(
        0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), EDist_CDF_Y_regrid.data());
    gitr_precision rfylAInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AphiDist_CDF_R_regrid.data());
    gitr_precision rfylAthetaInterpVal = interp3d(
        0.44, 5.0, std::log10(250.0), nA_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, angleDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), AthetaDist_CDF_R_regrid.data());
    gitr_precision rflEInterpVal = interp3d(
        0.44, 63.0, std::log10(10.0), nE_sputtRefDistOut, nA_sputtRefDistIn,
        nE_sputtRefDistIn, energyDistGrid01.data(), A_sputtRefDistIn.data(),
        Elog_sputtRefDistIn.data(), EDist_CDF_R_regrid.data());

    std::cout << "Finished surface model import sputtering" << spylAInterpVal
                << " " << spylAthetaInterpVal << " " << sputEInterpVal
                << std::endl;
    std::cout << "Finished surface model import reflection" << rfylAInterpVal
                << " " << rfylAthetaInterpVal << " " << rflEInterpVal
                << std::endl;
    
    // Put all the data into a tuple --> will eventually be a struct
    auto result = std::make_tuple(nE_sputtRefCoeff, nA_sputtRefCoeff, E_sputtRefCoeff, A_sputtRefCoeff, Elog_sputtRefCoeff, energyDistGrid01, angleDistGrid01, spyl_surfaceModel, rfyl_surfaceModel,
        nE_sputtRefDistIn, nA_sputtRefDistIn,
        E_sputtRefDistIn, A_sputtRefDistIn, Elog_sputtRefDistIn, E_sputtRefDistOut, E_sputtRefDistOutRef, Aphi_sputtRefDistOut, Atheta_sputtRefDistOut, EDist_Y, AphiDist_Y, AthetaDist_Y,
        EDist_R, AphiDist_R, AthetaDist_R, nDistE_surfaceModel, nDistA_surfaceModel, nDistE_surfaceModelRef, nE_sputtRefDistOut, nA_sputtRefDistOut,
        energyDistGrid01Ref, EDist_CDF_Y, AphiDist_CDF_Y, AthetaDist_CDF_Y, EDist_CDF_R, AphiDist_CDF_R, AthetaDist_CDF_R, EDist_CDF_Y_regrid, AphiDist_CDF_Y_regrid,
        AthetaDist_CDF_Y_regrid, EDist_CDF_R_regrid, AphiDist_CDF_R_regrid, AthetaDist_CDF_R_regrid, spylAInterpVal, spylAthetaInterpVal, sputEInterpVal, rfylAInterpVal, rfylAthetaInterpVal, rflEInterpVal,
        nEdist, nAdist, E0dist, Edist, A0dist, Adist, nE_sputtRefDistOutRef);
    
    cacheSurfaceData[key] = result;
    return result;
}

