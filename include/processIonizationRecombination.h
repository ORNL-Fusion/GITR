//------------------------------------------------------------------------------
// GITR: processIonizationRecombination.h
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - August 2023 by Diaw
//
// Description:
//     Processes ADAS data for ionization and recombination. Requires an ADAS
//     file for each material with a nuclear charge Z. The file should be named
//     in the format: ADAS_Rates_Z.nc and be located in the input/adasData 
//     directory.
// FIXME: we eventually wants to move all this data into a single file.
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------

#include <netcdf>

enum ElementaryProcess
{
    IONIZATION,
    RECOMBINATION
};


/* Function to process ADAS data for the elementary processes listed above */
inline std::tuple<size_t, size_t, size_t, std::vector<double>, std::vector<double>, std::vector<double>> process_rates(ElementaryProcess process, int charge)
{

    std::string filePrefix;
    if (process == IONIZATION)
    {
        filePrefix = "Ionization";
    }
    else if (process == RECOMBINATION)
    {
        filePrefix = "Recombination";
    }
    std::string input_path = "input/adasData/";
    std::string ratesFiles = "ADAS_Rates_" + std::to_string(charge) + ".nc";
    
    // Open the netCDF file
    netCDF::NcFile data(input_path + ratesFiles, netCDF::NcFile::read);


    // Get the variable objects
    netCDF::NcVar gridTemperature_var = data.getVar("gridTemperature_" + filePrefix);
    netCDF::NcVar gridDensity_var = data.getVar("gridDensity_" + filePrefix);
    netCDF::NcVar rateCoeff_var = data.getVar(filePrefix + "RateCoeff");

    // Get the dimensions of the variables
    std::vector<netCDF::NcDim> temperatureDims = gridTemperature_var.getDims();
    std::vector<netCDF::NcDim> densityDims = gridDensity_var.getDims();
    std::vector<netCDF::NcDim> chargeStateDims = rateCoeff_var.getDims();

    // Get the number of temperatures, densities, and charge states
    size_t nTemperatures = temperatureDims[0].getSize();
    size_t nDensities = densityDims[0].getSize();
    size_t nChargeStates = chargeStateDims[0].getSize();

    // Resize the arrays
    std::vector<double> gridTemperature(nTemperatures);
    std::vector<double> gridDensity(nDensities);
    std::vector<double> rateCoeff(nChargeStates * nTemperatures * nDensities);

    // Read the variable data into the arrays
    gridTemperature_var.getVar(gridTemperature.data());
    gridDensity_var.getVar(gridDensity.data());
    rateCoeff_var.getVar(rateCoeff.data());

    // Close the netCDF file
    data.close();

    return std::make_tuple(nTemperatures, nDensities, nChargeStates, gridTemperature, gridDensity, rateCoeff);
}
