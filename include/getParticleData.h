#include <iostream>
#include <fstream>
#include <string>
#include <netcdf.h>
#include <vector>
#include <libconfig.h++>
#include "utils.h"

const gitr_precision PI = 3.141592653589793;
const gitr_precision ELECTRON_CHARGE = 1.602e-19;
const gitr_precision AMU_CONVERSION = 1.66e-27;

template<typename T>
std::vector<T> getVarData(const netCDF::NcVar& var) {
    std::vector<T> data(var.getDim(0).getSize());
    var.getVar(&data[0]);
    return data;
}

std::vector<std::string> getVarData(netCDF::NcFile& file, const std::string& varName) {
    int varId;
    nc_inq_varid(file.getId(), varName.c_str(), &varId);
    nc_type varType;
    nc_inq_vartype(file.getId(), varId, &varType);
    if (varType != NC_STRING) {
        throw std::runtime_error("The variable is not of string type.");
    }
    size_t varLen;
    nc_inq_dimlen(file.getId(), 0, &varLen); // Assuming 0 is the dimension ID for 'nP'
    std::vector<char*> buffer(varLen);
    nc_get_var_string(file.getId(), varId, buffer.data());
    std::vector<std::string> result(varLen);
    for (size_t i = 0; i < varLen; ++i) {
        result[i] = buffer[i];
        free(buffer[i]); // Free memory allocated by nc_get_var_string
    }
    return result;
}


struct ParticleData {
    std::vector<gitr_precision> xpfile, ypfile, zpfile,
        vxpfile, vypfile, vzpfile, charge, amu, Z, species;
};


ParticleData readParticleData(const std::string &ncParticleSourceFile) {
    // Read in particle source file
    ParticleData data;

    std::cout << "Starting particleSourcefile import " << std::endl;

    netCDF::NcFile ncp;
    try {
        ncp.open("input/" + ncParticleSourceFile, netCDF::NcFile::read);
    } catch (netCDF::exceptions::NcException &e) {
        std::cerr << e.what() << std::endl;
        std::cout << "FAILURE*************************************" << std::endl;
        throw std::runtime_error("Failed to open NetCDF file.");
    }

    try {
        // Fetch the dimension
        netCDF::NcDim ps_nP = ncp.getDim("nP");
        size_t nPfile = ps_nP.getSize();

        // Resize the data vectors
        data.xpfile.resize(nPfile); 
        data.ypfile.resize(nPfile);
        data.zpfile.resize(nPfile);
        data.vxpfile.resize(nPfile);
        data.vypfile.resize(nPfile);
        data.vzpfile.resize(nPfile);
        data.charge.resize(nPfile);
        data.amu.resize(nPfile);
        data.Z.resize(nPfile);
        data.species.resize(nPfile);

        // Fetch the variables
        data.xpfile = getVarData<gitr_precision>(ncp.getVar("x"));
        data.ypfile = getVarData<gitr_precision>(ncp.getVar("y"));
        data.zpfile = getVarData<gitr_precision>(ncp.getVar("z"));
        data.vxpfile = getVarData<gitr_precision>(ncp.getVar("vx"));
        data.vypfile = getVarData<gitr_precision>(ncp.getVar("vy"));
        data.vzpfile = getVarData<gitr_precision>(ncp.getVar("vz"));
        data.charge = getVarData<gitr_precision>(ncp.getVar("charge"));
        data.amu = getVarData<gitr_precision>(ncp.getVar("amu"));
        data.Z = getVarData<gitr_precision>(ncp.getVar("Z"));
        data.species = getVarData<gitr_precision>(ncp.getVar("species"));

    } catch (netCDF::exceptions::NcException &e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("Error while reading data from NetCDF file.");
    }

    ncp.close();
    return data;
    }


ParticleData generateParticleData(const libconfig::Config &cfg, int particle_source_energy) {
    ParticleData data;

    std::cout << "creating particle source" << std::endl;
    int nP = cfg.lookup("impurityParticleSource.nP");

    data.xpfile.resize(nP); data.ypfile.resize(nP); data.zpfile.resize(nP);
    data.vxpfile.resize(nP); data.vypfile.resize(nP); data.vzpfile.resize(nP);
    data.charge.resize(nP); data.amu.resize(nP); data.Z.resize(nP);
    data.species.resize(nP);

    libconfig::Setting& speciesArray = cfg.lookup("impurityParticleSource.initialConditions.species");
    int num_species = speciesArray.getLength();
    printf("Number of species: %d\n", num_species);
    std::vector<gitr_precision> amu(num_species), Z(num_species), charge(num_species);
    int particlesPerSpecies[num_species];
    int particle_index = 0;

    for (int i = 0; i < num_species; i++) {
        gitr_precision _amu, _Z, _charge;
                 gitr_precision _energy_eV, _phi, _theta;
                 gitr_precision _x_start, _y_start, _z_start;
        
                 int  _numberParticlesPerSpecies;
                 if (speciesArray[i].lookupValue("impurity_amu", _amu) &&
                     speciesArray[i].lookupValue("impurity_Z", _Z) &&
                     speciesArray[i].lookupValue("charge", _charge) &&
                     speciesArray[i].lookupValue("numberParticlesPerSpecies", _numberParticlesPerSpecies) &&
                     speciesArray[i].lookupValue("energy_eV", _energy_eV) &&
                     speciesArray[i].lookupValue("theta", _theta) &&
                     speciesArray[i].lookupValue("x_start", _x_start) &&
                     speciesArray[i].lookupValue("y_start", _y_start) &&
                     speciesArray[i].lookupValue("z_start", _z_start) &&
                     speciesArray[i].lookupValue("phi", _phi) )
                      {
                     // Store the values
                     amu[i] = _amu;
                     Z[i] = _Z;
                     charge[i] = _charge;
                     particlesPerSpecies[i] = _numberParticlesPerSpecies;
                     std::cout << "Species " << (i+1) << " - Impurity amu Z charge: "
                               << _amu << " " << _Z << " " << _charge << " particlesPerSpecies"  << " " << _numberParticlesPerSpecies << " "
                               << " phi theta energy: " << _phi << " " << _theta << " " << _energy_eV << std::endl;
                 } else {
                     std::cout << "ERROR: Could not get initial conditions for species "
                               << (i+1) << std::endl;
                 }
               
        if (particle_source_energy == 0) {
            for (int j = 0; j < _numberParticlesPerSpecies; j++) {
                gitr_precision rad_phi = _phi * M_PI / 180.0;
                gitr_precision rad_theta = _theta * M_PI / 180.0;
                
                gitr_precision vtotal = std::sqrt(2.0 * _energy_eV * ELECTRON_CHARGE/ _amu / AMU_CONVERSION);
                
                gitr_precision vx = vtotal * std::sin(rad_phi) * std::cos(rad_theta);
                gitr_precision vy = vtotal * std::sin(rad_phi) * std::sin(rad_theta);
                gitr_precision vz = vtotal * std::cos(rad_phi);
                
                if (rad_phi == 0.0) {
                    vx = 0.0;
                    vy = 0.0;
                    vz = vtotal;
                }
                data.xpfile[particle_index] = _x_start;
                data.ypfile[particle_index] = _y_start;
                data.zpfile[particle_index] = _z_start;
                data.vxpfile[particle_index] = vx;
                data.vypfile[particle_index] = vy;
                data.vzpfile[particle_index] = vz;
                data.charge[particle_index] = _charge;
                data.amu[particle_index] = _amu;
                data.Z[particle_index] = _Z;
                data.species[particle_index] = i;
                
                particle_index++;
            }
        } else {
                        std::cout << "Warning: Non-zero particle_source_energy is not currently supported." << std::endl;
        }
    }
    std::cout << "size of data: " << data.xpfile.size() << std::endl;
    return data;
}

 void initializeParticleArray(const ParticleData& particleData, Particles* particleArray,
    sim::Array<gitr_precision>& px, sim::Array<gitr_precision>& py, sim::Array<gitr_precision>& pz,
    sim::Array<gitr_precision>& pvx, sim::Array<gitr_precision>& pvy, sim::Array<gitr_precision>& pvz,
    sim::Array<gitr_precision>& pZ, sim::Array<gitr_precision>& pamu, sim::Array<gitr_precision>& pcharge,
    gitr_precision dt, sim::Array<int >& pSpecies)
    {

    auto& xpfile = particleData.xpfile;
    auto& ypfile = particleData.ypfile;
    auto& zpfile = particleData.zpfile;
    auto& vxpfile = particleData.vxpfile;
    auto& vypfile = particleData.vypfile;
    auto& vzpfile = particleData.vzpfile;
    auto& charge = particleData.charge;
    auto& mass = particleData.amu;
    auto& ionizationState = particleData.Z;
    auto& _species = particleData.species;

    long nParticles = particleData.Z.size();

    for (size_t i = 0; i < nParticles; ++i) {
        // Set particle data
        gitr_precision x = xpfile[i];
        gitr_precision y = ypfile[i];
        gitr_precision z = zpfile[i];
        gitr_precision vx = vxpfile[i];
        gitr_precision vy = vypfile[i];
        gitr_precision vz = vzpfile[i];
        gitr_precision particleCharge = charge[i];
        gitr_precision particleMass = mass[i];
        gitr_precision ionization = ionizationState[i];
        int species = _species[i];

        particleArray->setParticleV(i, x, y, z, vx, vy, vz, ionization, particleMass, particleCharge, dt, species);
 
        px[i] = x;
        py[i] = y;
        pz[i] = z;
        pvx[i] = vx;
        pvy[i] = vy;
        pvz[i] = vz;
        pZ[i] = ionization;
        pamu[i] = particleMass;
        pcharge[i] = particleCharge;
        pSpecies[i] = species;

    }

}