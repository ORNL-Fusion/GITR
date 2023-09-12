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
//
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
    int nSurfaces, Surfaces* surfaces, int nEdist, int nAdist)
{
// add initial particle erosion to surface counting

        int nSpecies = 0;
        std::vector<int> surfaceNumbers(nSurfaces, 0);
        int srf = 0;
        for (int i = 0; i < nLines; i++) {
        if (boundaries[i].surface) {
            surfaceNumbers[srf] = i;

            surfaces->grossErosion[srf] = surfaces->grossErosion[srf] + grossErosion[srf];
            srf = srf + 1;
        }
        }  
        netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);
        netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
        vector<netCDF::NcDim> dims1;
        dims1.push_back(nc_nLines);

        vector<netCDF::NcDim> dimsSurfE;
        dimsSurfE.push_back(nc_nLines);
        netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
        netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
        dimsSurfE.push_back(nc_nEnergies);
        dimsSurfE.push_back(nc_nAngles);
        netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, nc_nLines);
        netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, nc_nLines);
        netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, nc_nLines);
        netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, nc_nLines);
        netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, nc_nLines);
        netCDF::NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, nc_nLines);
        netCDF::NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike", netcdf_precision, nc_nLines);
        nc_grossDep.putVar(&surfaces->grossDeposition[0]);
        nc_surfNum.putVar(&surfaceNumbers[0]);
        nc_grossEro.putVar(&surfaces->grossErosion[0]);
        nc_aveSpyl.putVar(&surfaces->aveSputtYld[0]);
        nc_spylCounts.putVar(&surfaces->sputtYldCount[0]);
        nc_sumParticlesStrike.putVar(&surfaces->sumParticlesStrike[0]);
        nc_sumWeightStrike.putVar(&surfaces->sumWeightStrike[0]);
        netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
        netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
        netCDF::NcVar nc_surfSputtDist = ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);
        nc_surfEDist.putVar(&surfaces->energyDistribution[0]);
        nc_surfReflDist.putVar(&surfaces->reflDistribution[0]);
        nc_surfSputtDist.putVar(&surfaces->sputtDistribution[0]);
        ncFile1.close();
    }


// void writeSurfaceData( int surface_model, 
//         {
//     if( surface_model > 0 || flux_ea > 0 )
//     {
//         int nSpecies = 0;
//         std::vector<int> surfaceNumbers(nSurfaces, 0);
//         int srf = 0;
//         for (int i = 0; i < nLines; i++) {
//         if (boundaries[i].surface) {
//             surfaceNumbers[srf] = i;

//             surfaces->grossErosion[srf] = surfaces->grossErosion[srf] + grossErosion[srf];
//             srf = srf + 1;
//         }
//         }  
//         netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);
//         netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
//         vector<netCDF::NcDim> dims1;
//         dims1.push_back(nc_nLines);

//         vector<netCDF::NcDim> dimsSurfE;
//         dimsSurfE.push_back(nc_nLines);
//         netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
//         netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
//         dimsSurfE.push_back(nc_nEnergies);
//         dimsSurfE.push_back(nc_nAngles);
//         netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, nc_nLines);
//         netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, nc_nLines);
//         netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, nc_nLines);
//         netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, nc_nLines);
//         netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, nc_nLines);
//         netCDF::NcVar nc_sumParticlesStrike = ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, nc_nLines);
//         netCDF::NcVar nc_sumWeightStrike = ncFile1.addVar("sumWeightStrike", netcdf_precision, nc_nLines);
//         nc_grossDep.putVar(&surfaces->grossDeposition[0]);
//         nc_surfNum.putVar(&surfaceNumbers[0]);
//         nc_grossEro.putVar(&surfaces->grossErosion[0]);
//         nc_aveSpyl.putVar(&surfaces->aveSputtYld[0]);
//         nc_spylCounts.putVar(&surfaces->sputtYldCount[0]);
//         nc_sumParticlesStrike.putVar(&surfaces->sumParticlesStrike[0]);
//         nc_sumWeightStrike.putVar(&surfaces->sumWeightStrike[0]);
//         netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
//         netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
//         netCDF::NcVar nc_surfSputtDist = ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);
//         nc_surfEDist.putVar(&surfaces->energyDistribution[0]);
//         nc_surfReflDist.putVar(&surfaces->reflDistribution[0]);
//         nc_surfSputtDist.putVar(&surfaces->sputtDistribution[0]);
//         ncFile1.close();
//     }
// }
