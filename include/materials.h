#ifndef MATERIALS_H
#define MATERIALS_H

#include <map>
#include <string>

#if USE_DOUBLE
using gitr_precision = double;
#else
using gitr_precision = float;
#endif

struct MaterialProperties {
    std::string name;
    gitr_precision surfaceBindingEnergy;
    gitr_precision mass;

    MaterialProperties() : name(""), surfaceBindingEnergy(0.0), mass(0.0) {}

    MaterialProperties(const std::string& name_, gitr_precision surfaceBindingEnergy_, gitr_precision mass_)
        : name(name_), surfaceBindingEnergy(surfaceBindingEnergy_), mass(mass_) {}
};

// Extern declarations
extern std::map<int, MaterialProperties> materialData;
extern std::map<gitr_precision, MaterialProperties> materialDataByMass;

bool materialExists(int nuclearCharge);
bool materialExistsByMass(gitr_precision mass);

#endif // MATERIALS_H
