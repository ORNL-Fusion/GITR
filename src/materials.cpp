#include "materials.h"

// Define the variable for nuclear charge
std::map<int, MaterialProperties> materialData = {
    {1,   MaterialProperties("H", 0.1, 1.008)},
    {2,   MaterialProperties("He", 0.2, 4.0026)},
    {3,   MaterialProperties("Li", 0.3, 6.94)},
    {6,   MaterialProperties("C", 0.6, 12.011)},
    {7,   MaterialProperties("N", 0.7, 14.007)},
    {8,   MaterialProperties("O", 0.8, 15.999)},
    {10,  MaterialProperties("Ne", 1.0, 20.180)},
    {12,  MaterialProperties("Mg", 1.2, 24.305)},
    {14,  MaterialProperties("Si", 1.4, 28.085)},
    {16,  MaterialProperties("S", 1.6, 32.06)},
    {18,  MaterialProperties("Ar", 1.8, 39.948)},
    {20,  MaterialProperties("Ca", 2.0, 40.078)},
    {22,  MaterialProperties("Ti", 2.2, 47.867)},
    {24,  MaterialProperties("Cr", 2.4, 51.996)},
    {26,  MaterialProperties("Fe", 2.6, 55.845)},
    {28,  MaterialProperties("Ni", 2.8, 58.693)},
    {30,  MaterialProperties("Zn", 3.0, 65.38)},
    {32,  MaterialProperties("Ge", 3.2, 72.63)},
    {34,  MaterialProperties("Se", 3.4, 78.971)},
    {36,  MaterialProperties("Kr", 3.6, 83.798)},
    {38,  MaterialProperties("Sr", 3.8, 87.62)},
    {40,  MaterialProperties("Zr", 4.0, 91.224)},
    {42,  MaterialProperties("Mo", 4.2, 95.95)},
    {44,  MaterialProperties("Ru", 4.4, 101.07)},
    {46,  MaterialProperties("Pd", 4.6, 106.42)},
    {48,  MaterialProperties("Cd", 4.8, 112.41)},
    {50,  MaterialProperties("Sn", 5.0, 118.71)},
    {52,  MaterialProperties("Te", 5.2, 127.60)},
    {74, MaterialProperties("W", 7.4, 183.84)},
    {76,  MaterialProperties("Os", 7.6, 190.23)},
    {78,  MaterialProperties("Pt", 7.8, 195.08)},
    {80,  MaterialProperties("Hg", 8.0, 200.59)},
    {82,  MaterialProperties("Pb", 8.2, 207.2)},
    {92,  MaterialProperties("U", 9.2, 238.03)},
    {94,  MaterialProperties("Pu", 9.4, 244.06)},
    {96,  MaterialProperties("Cm", 9.6, 247.07)},
    {98,  MaterialProperties("Cf", 9.8, 251.08)},
    {100, MaterialProperties("Fm", 10.0, 252.08)},
    {102, MaterialProperties("No", 10.2, 257.10)},
    {104, MaterialProperties("Rf", 10.4, 261.11)},
    {106, MaterialProperties("Sg", 10.6, 262.11)},
    {108, MaterialProperties("Hs", 10.8, 265.12)},
    {110, MaterialProperties("Ds", 11.0, 268.13)},
    {112, MaterialProperties("Cn", 11.2, 271.13)},
    {114, MaterialProperties("Fl", 11.4, 274.14)},
    {116, MaterialProperties("Lv", 11.6, 277.15)},
    {118, MaterialProperties("Og", 11.8, 281.16)}
};


// Define the variable for mass
std::map<gitr_precision, MaterialProperties> materialDataByMass;

void initializeMaterialsByMass() {
    for (const auto& pair : materialData) {
        materialDataByMass[pair.second.mass] = pair.second;
    }
}

bool materialExists(int nuclearCharge) {
    return materialData.find(nuclearCharge) != materialData.end();
}

bool materialExistsByMass(gitr_precision mass) {
    return materialDataByMass.find(mass) != materialDataByMass.end();
}