//------------------------------------------------------------------------------
// GITR: diagnostics.cpp
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - September 2023 Diaw
//
// Description:
//     This source file handles the storage of particle data in a netCDF format.
//     It interfaces with the netCDF library and ensures the data structure and 
//     variable names conform to the desired GITR data output standards.
//
// Functions Included:
//     1. storeParticleData:
//        - Takes in a base filename, particle array, and the number of particles.
//        - Writes particle data to a netCDF file format.
//        - Variables include particle position, velocity, transit time, and more.
//     2. ParticleErosionToSurface:
//        - Takes in a boundary array, number of boundaries, gross erosion array,
//          number of surfaces, surface array, number of energy bins, and number
//          of angle bins.
//        - Writes surface data to a netCDF file format.
//        - Variables include gross deposition, gross erosion, average sputter
//          yield, and more.
//     3. writeParticleDataHistories:
//        - Takes in particle position, velocity, charge, atomic number, weight,
//          number of histories per particle, and number of particles.
//        - Writes particle history data to a netCDF file format.
//        - Variables include particle position, velocity, charge, atomic number,
//          weight, and more.
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------

#include "Particles.h"
#include "array.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <netcdf>
#include <Surfaces.h>
#include "Boundary.h"

extern netCDF::NcType netcdf_precision;

void storeParticleData(const std::string& base_filename, Particles *particleArray, int nP){

    // Write netCDF output for positions
    std::string filename = "output/" + base_filename;

   // Create and set dimensions for netCDF output
    netCDF::NcFile ncFile0(filename, netCDF::NcFile::replace);
    netCDF::NcDim nc_nP0 = ncFile0.addDim("nP", nP);
    std::vector<netCDF::NcDim> dims0 = { nc_nP0 };

    netCDF::NcVar nc_x0 = ncFile0.addVar("x", netcdf_precision, dims0);
    netCDF::NcVar nc_y0 = ncFile0.addVar("y", netcdf_precision, dims0);
    netCDF::NcVar nc_z0 = ncFile0.addVar("z", netcdf_precision, dims0);
    netCDF::NcVar nc_vx0 = ncFile0.addVar("vx", netcdf_precision, dims0);
    netCDF::NcVar nc_vy0 = ncFile0.addVar("vy", netcdf_precision, dims0);
    netCDF::NcVar nc_vz0 = ncFile0.addVar("vz", netcdf_precision, dims0);
    netCDF::NcVar nc_trans0 = ncFile0.addVar("transitTime", netcdf_precision, dims0);
    netCDF::NcVar nc_impact0 = ncFile0.addVar("hitWall", netcdf_precision, dims0);
    netCDF::NcVar nc_surfHit0 = ncFile0.addVar("surfaceHit", netCDF::ncInt, dims0);
    netCDF::NcVar nc_weight0 = ncFile0.addVar("weight", netcdf_precision, dims0);
    netCDF::NcVar nc_charge0 = ncFile0.addVar("charge", netcdf_precision, dims0);
    netCDF::NcVar nc_leak0 = ncFile0.addVar("hasLeaked", netCDF::ncInt, dims0);
    netCDF::NcVar nc_dist0 = ncFile0.addVar("distTraveled", netcdf_precision, dims0);
    netCDF::NcVar nc_time0 = ncFile0.addVar("time", netcdf_precision, dims0);
    netCDF::NcVar nc_dt0 = ncFile0.addVar("dt", netcdf_precision, dims0);
    netCDF::NcVar nc_mass0 = ncFile0.addVar("amu", netcdf_precision, dims0);
    netCDF::NcVar nc_Z0 = ncFile0.addVar("Z", netcdf_precision, dims0);
    netCDF::NcVar nc_species0 = ncFile0.addVar("species", netCDF::ncInt, dims0);

    nc_x0.putVar(&particleArray->xprevious[0]);
    nc_y0.putVar(&particleArray->yprevious[0]);
    nc_z0.putVar(&particleArray->zprevious[0]);
    nc_vx0.putVar(&particleArray->vx[0]);
    nc_vy0.putVar(&particleArray->vy[0]);
    nc_vz0.putVar(&particleArray->vz[0]);
    nc_trans0.putVar(&particleArray->transitTime[0]);
    nc_impact0.putVar(&particleArray->hitWall[0]);
    nc_surfHit0.putVar(&particleArray->surfaceHit[0]);
    nc_weight0.putVar(&particleArray->weight[0]);
    nc_charge0.putVar(&particleArray->charge[0]);
    nc_leak0.putVar(&particleArray->hasLeaked[0]);
    nc_dist0.putVar(&particleArray->distTraveled[0]);
    nc_time0.putVar(&particleArray->time[0]);
    nc_dt0.putVar(&particleArray->dt[0]);
    nc_mass0.putVar(&particleArray->amu[0]);
    nc_Z0.putVar(&particleArray->Z[0]);
    nc_species0.putVar(&particleArray->species[0]);
    ncFile0.close();

}


void ParticleErosionToSurface(Boundary* boundaries, int nLines, gitr_precision* grossErosion, 
    int nSurfaces, Surfaces* surfaces, int nEdist, int nAdist, int nSpecies)
{
    // Create a vector for surface numbers
    std::vector<int> surfaceNumbers(nSurfaces, 0);
    int srf = 0;
    for (int i = 0; i < nLines; i++) {
        if (boundaries[i].surface) {
                surfaceNumbers[srf] = i;
                surfaces->grossErosion[srf] = surfaces->grossErosion[srf] + grossErosion[srf];

                srf++;
            }

    }  
    try {
        netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);

        // Define dimensions
        netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
        std::vector<netCDF::NcDim> dims1 = {nc_nLines};

        // Create more dimensions for energy and angle distribution
        netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
        netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
        std::vector<netCDF::NcDim> dimsSurfE = {nc_nLines, nc_nEnergies, nc_nAngles};

        // Define variables
        netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, dims1);
        netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, dims1);
        netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, dims1);
        netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, dims1);
        netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, dims1);
        netCDF::NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, dims1);
        netCDF::NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike", netcdf_precision, dims1);
        netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
        netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
        netCDF::NcVar nc_surfSputtDist = ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);

        // Write data to file
        nc_grossDep.putVar(&surfaces->grossDeposition[0]);
        nc_surfNum.putVar(surfaceNumbers.data());
        nc_grossEro.putVar(&surfaces->grossErosion[0]);
        nc_aveSpyl.putVar(&surfaces->aveSputtYld[0]);
        nc_spylCounts.putVar(&surfaces->sputtYldCount[0]);
        nc_sumParticlesStrike.putVar(&surfaces->sumParticlesStrike[0]);
        nc_sumWeightStrike.putVar(&surfaces->sumWeightStrike[0]);
        nc_surfEDist.putVar(&surfaces->energyDistribution[0]);
        nc_surfReflDist.putVar(&surfaces->reflDistribution[0]);
        nc_surfSputtDist.putVar(&surfaces->sputtDistribution[0]);

        // Close the NetCDF file
        ncFile1.close();
    }
    catch (netCDF::exceptions::NcException& e) {
        std::cerr << "NetCDF Error: " << e.what() << std::endl;
    }
}



void writeParticleDataHistories(
        const sim::Array<gitr_precision>& positionHistoryX, 
        const sim::Array<gitr_precision>& positionHistoryY, 
        const sim::Array<gitr_precision>& positionHistoryZ,
        const sim::Array<gitr_precision>& velocityHistoryX, 
        const sim::Array<gitr_precision>& velocityHistoryY, 
        const sim::Array<gitr_precision>& velocityHistoryZ, 
        const sim::Array<gitr_precision>& chargeHistory,
        const sim::Array<gitr_precision>& ZHistory,
        const sim::Array<gitr_precision>& weightHistory,
        const sim::Array<gitr_precision>& massHistory,
        const sim::Array<gitr_precision>& surfaceHitHistory,
        const sim::Array<gitr_precision>& hitWallHistory,
        int nHistoriesPerParticle, 
        int nP
    ) {
        // Write netCDF output for histories
        netCDF::NcFile ncFile_hist("output/history.nc", netCDF::NcFile::replace);

        netCDF::NcDim nc_nT = ncFile_hist.addDim("nT", nHistoriesPerParticle);
        netCDF::NcDim nc_nP = ncFile_hist.addDim("nP", nP);
        
        std::vector<netCDF::NcDim> dims_hist;
        dims_hist.push_back(nc_nP);
        dims_hist.push_back(nc_nT);

        netCDF::NcVar nc_x = ncFile_hist.addVar("x", netcdf_precision, dims_hist);
        netCDF::NcVar nc_y = ncFile_hist.addVar("y", netcdf_precision, dims_hist);
        netCDF::NcVar nc_z = ncFile_hist.addVar("z", netcdf_precision, dims_hist);

        netCDF::NcVar nc_vx = ncFile_hist.addVar("vx", netcdf_precision, dims_hist);
        netCDF::NcVar nc_vy = ncFile_hist.addVar("vy", netcdf_precision, dims_hist);
        netCDF::NcVar nc_vz = ncFile_hist.addVar("vz", netcdf_precision, dims_hist);

        netCDF::NcVar nc_charge = ncFile_hist.addVar("charge", netcdf_precision, dims_hist);
        netCDF::NcVar nc_Z = ncFile_hist.addVar("Z", netcdf_precision, dims_hist);
        netCDF::NcVar nc_weight = ncFile_hist.addVar("weight", netcdf_precision, dims_hist);
        netCDF::NcVar nc_mass = ncFile_hist.addVar("amu", netcdf_precision, dims_hist);
        netCDF::NcVar nc_surfaceHit = ncFile_hist.addVar("surfaceHit", netCDF::ncInt, dims_hist);
        netCDF::NcVar nc_hitWall = ncFile_hist.addVar("hitWall", netCDF::ncInt, dims_hist);


        nc_x.putVar(&positionHistoryX.data()[0]);
        nc_y.putVar(&positionHistoryY.data()[0]);
        nc_z.putVar(&positionHistoryZ.data()[0]);

        nc_vx.putVar(&velocityHistoryX.data()[0]);
        nc_vy.putVar(&velocityHistoryY.data()[0]);
        nc_vz.putVar(&velocityHistoryZ.data()[0]);

        nc_charge.putVar(&chargeHistory.data()[0]);
        nc_Z.putVar(&ZHistory.data()[0]);
        nc_weight.putVar(&weightHistory.data()[0]);
        nc_mass.putVar(&massHistory.data()[0]);

        nc_surfaceHit.putVar(&surfaceHitHistory.data()[0]);
        nc_hitWall.putVar(&hitWallHistory.data()[0]);

        ncFile_hist.close();
    }

void writeSpecData(
    int nBins, 
    int net_nX, 
    int net_nY, 
    int net_nZ, 
    int spectroscopy, 
    const sim::Array<gitr_precision>&gridX_bins, 
    const sim::Array<gitr_precision>&gridY_bins, 
    const sim::Array<gitr_precision>&gridZ_bins, 
    const sim::Array<gitr_precision>&net_Bins
) {
    // Write netCDF output for density data
    netCDF::NcFile ncFile("output/spec.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nBins = ncFile.addDim("nBins", nBins + 1);
    netCDF::NcDim nc_nR = ncFile.addDim("nR", net_nX);
    netCDF::NcDim nc_nY;
    if(spectroscopy > 2) {
        nc_nY = ncFile.addDim("nY", net_nY);
    }

    netCDF::NcDim nc_nZ = ncFile.addDim("nZ", net_nZ);
    std::vector<netCDF::NcDim> dims;
    dims.push_back(nc_nBins);
    dims.push_back(nc_nZ);

    if(spectroscopy > 2) {
        dims.push_back(nc_nY);
    }

    dims.push_back(nc_nR);
    netCDF::NcVar nc_n = ncFile.addVar("n", netcdf_precision, dims);
    netCDF::NcVar nc_gridR = ncFile.addVar("gridR", netcdf_precision, nc_nR); 
    netCDF::NcVar nc_gridZ = ncFile.addVar("gridZ", netcdf_precision, nc_nZ); 

    nc_gridR.putVar(&gridX_bins[0]);
    nc_gridZ.putVar(&gridZ_bins[0]);

    if(spectroscopy > 2) {
        netCDF::NcVar nc_gridY = ncFile.addVar("gridY", netcdf_precision, nc_nY); 
        nc_gridY.putVar(&gridY_bins[0]);
    }
    
    nc_n.putVar(&net_Bins[0]);
    ncFile.close();
}

void writeParticleSourceFile(
    int nP,
   const sim::Array<gitr_precision>& pvx,
   const sim::Array<gitr_precision>& pvy,
   const sim::Array<gitr_precision>& pvz,
   const sim::Array<gitr_precision>& px,
   const sim::Array<gitr_precision>& py,
   const sim::Array<gitr_precision>& pz,
   const sim::Array<gitr_precision>& pcharge,
   const sim::Array<gitr_precision>& pamu,
   const sim::Array<gitr_precision>& pZ,
   const sim::Array<int>& speciesType
) {
    std::cout << "writing particles out file" << std::endl;

    netCDF::NcFile ncFile_particles("output/particleSource.nc", netCDF::NcFile::replace);
    netCDF::NcDim pNP = ncFile_particles.addDim("nP", nP);
    netCDF::NcVar p_vx = ncFile_particles.addVar("vx", netcdf_precision, pNP);
    netCDF::NcVar p_vy = ncFile_particles.addVar("vy", netcdf_precision, pNP);
    netCDF::NcVar p_vz = ncFile_particles.addVar("vz", netcdf_precision, pNP);
    netCDF::NcVar p_x = ncFile_particles.addVar("x", netcdf_precision, pNP);
    netCDF::NcVar p_y = ncFile_particles.addVar("y", netcdf_precision, pNP);
    netCDF::NcVar p_z = ncFile_particles.addVar("z", netcdf_precision, pNP);
    netCDF::NcVar p_charge = ncFile_particles.addVar("charge", netcdf_precision, pNP);
    netCDF::NcVar p_amu = ncFile_particles.addVar("amu", netcdf_precision, pNP);
    netCDF::NcVar p_Z = ncFile_particles.addVar("Z", netcdf_precision, pNP);
    netCDF::NcVar p_speciesType = ncFile_particles.addVar("speciesType", netCDF::ncInt, pNP);

    p_vx.putVar(&pvx[0]);
    p_vy.putVar(&pvy[0]);
    p_vz.putVar(&pvz[0]);
    p_x.putVar(&px[0]);
    p_y.putVar(&py[0]);
    p_z.putVar(&pz[0]);
    p_charge.putVar(&pcharge[0]);
    p_amu.putVar(&pamu[0]);
    p_Z.putVar(&pZ[0]);
    p_speciesType.putVar(&speciesType[0]);

    ncFile_particles.close();
}
