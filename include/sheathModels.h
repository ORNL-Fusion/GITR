//------------------------------------------------------------------------------
// GITR: sheathModels.h
//------------------------------------------------------------------------------
//
// Contributors:
//     - GITR Community
//
// Last Modified:
//     - August 2023 by Diaw
//
// Description:
//     This header file defines the models for sheath electric fields.
//
// Models Included:
//     1. Stangeby Model:
//        - Reference: 
//            - P.C. Stangeby "The Chodura sheath for angles of a few degrees 
//              between the magnetic field and the surface of divertor targets 
//              and limiters", Nucl. Fusion 52 083012 (2012).
//              [DOI: 10.1088/0029-5515/52/8/083012]
//              
//     2. Coulette-Manfredi Model:
//        - Reference:
//            - D. Coulette and G. Manfredi "Kinetic simulations of the Chodura 
//              and Debye sheaths for magnetic fields with grazing incidence", 
//              Plasma Phys. Control. Fusion 58 025008 (2016).
//              [DOI: 10.1088/0741-3335/58/2/025008]
//
// Note:
//     This file is a component of the GITR codebase.
//
//------------------------------------------------------------------------------


#include <vector>
#include <cmath>
#include "Boundary.h"


#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif


enum class SheathModel {
    STANGEBY = 0,
    COULETTE_MANFREDI = 1
};

class ISheathModel {
public:
    virtual ~ISheathModel() {}  
    virtual gitr_precision getElectricField(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) = 0;
};

class StangebyModel : public ISheathModel {
public:
    gitr_precision getElectricField(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) override {
        // Extract boundaries parameters 
        gitr_precision angle = boundaryVector[minIndex].angle;   
        gitr_precision fd  = boundaryVector[minIndex].fd;
        gitr_precision potential = boundaryVector[minIndex].potential;
        gitr_precision debyeLength = boundaryVector[minIndex].debyeLength;
        gitr_precision larmorRadius = boundaryVector[minIndex].larmorRadius;

        gitr_precision E_component_DS = fd / (2.0 * debyeLength) * std::exp(-minDistance / (2.0 * debyeLength));
        gitr_precision E_component_CD = (1.0 - fd) / larmorRadius * std::exp(-minDistance / larmorRadius);
        gitr_precision Emag = potential * (E_component_DS + E_component_CD);

        return Emag;
    }
};

class CouletteManfrediModel : public ISheathModel {
public:
    gitr_precision getElectricField(gitr_precision minDistance, Boundary* boundaryVector, int minIndex) override {
        // Extract boundaries parameters 
        gitr_precision alpha = boundaryVector[minIndex].angle; 
        
        gitr_precision C1_alpha = -0.00281407f * alpha - 2.31655435f;
        gitr_precision C2_alpha = 0.00640402f * alpha + 0.01023915f;

        gitr_precision te = boundaryVector[minIndex].te;
        gitr_precision debyeLength = boundaryVector[minIndex].debyeLength;

        gitr_precision thresholdForSheathModel =  200; // 200 debye lengths

        if ( minDistance/debyeLength > thresholdForSheathModel ){
            return 0;
        }
        else{
            gitr_precision result = C1_alpha * C2_alpha *  std::exp(-C2_alpha * abs(minDistance  / debyeLength )) * te / debyeLength ;
            return result;
        }
    }
};

gitr_precision sheathModel(int sheath_model_type, gitr_precision minDistance, Boundary* boundaryVector, int minIndex) {
    ISheathModel* sheath;
    
    switch (sheath_model_type) {
        case 0:
            sheath = new StangebyModel();
            break;
        case 1:
            sheath = new CouletteManfrediModel();
            break;
    }

    gitr_precision result = sheath->getElectricField(minDistance, boundaryVector, minIndex);
    delete sheath;  
    return result;
}

//END SHEATH MODEL