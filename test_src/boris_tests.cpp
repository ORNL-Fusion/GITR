#include <iostream>
#include <thrust/execution_policy.h>
#include "test_utils.hpp"
#include "config_interface.h"
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "boris.h"
#include "Surfaces.h"
#include "geometryCheck.h"
#include "boundaryInit.h"

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

TEST_CASE( "boris - not fully implemented" )
{
  SECTION( "getE tests" )
  {
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig( cfg_geom, GET_E_TEST_FILE );
    int nLines = 1;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    sim::Array<gitr_precision> ne(1, 1.0e19);
    
    // Temperature = 20 eV
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);
    
    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    /* required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1 */
    /* create a unified setup script */
    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    /* uniform bfield */
    br[ 0 ] = std::cos(M_PI*5.0/180);
    /* large bfield in teslas gives smaller gyromotion radius */
    by[ 0 ] = 0;
    bz[ 0 ] = -std::sin(M_PI*5.0/180);;

    /* for the uniform efield, set efield to 1000 in z just make the cross product geometry */
    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    gitr_precision background_Z = 1;
    gitr_precision background_amu = 2;
    gitr_precision biasPotential = 0;

  
    std::for_each(boundaries.begin(), boundaries.end() - 1,
                boundary_init(background_Z, background_amu, nR_Dens, nZ_Dens,
                              DensGridr.data(), DensGridz.data(), ni.data(),
                              ne.data(), nR_Bfield, nZ_Bfield,
                              bfieldGridr.data(), bfieldGridz.data(), br.data(),
                              bz.data(), by.data(), nR_Temp, nZ_Temp,
                              TempGridr.data(), TempGridz.data(), ti.data(),
                              te.data(), biasPotential));
    
    int nHashes = 1;
    int nR_closeGeom_sheath = 1;
    int nY_closeGeom_sheath = 1;
    int nZ_closeGeom_sheath = 1;
    int nHashPoints_sheath = 1;
    int n_closeGeomElements_sheath = 1;
    sim::Array<gitr_precision> closeGeomGridr_sheath(1),
      closeGeomGridy_sheath(1), closeGeomGridz_sheath(1);
    sim::Array<int> closeGeom_sheath(1, 0);
    
    int closestBoundaryIndex = 0;
    int surfIndex = 0;
    gitr_precision minDistance = 0.0;
    gitr_precision thisE[3] = {0.0};
    sim::Array<gitr_precision> px(1, 0);
    sim::Array<gitr_precision> py(1, 0);
    sim::Array<gitr_precision> pz(1, 0.001);
    gitr_precision dz = 0.005/10000.0;
     int nZ = 10000;
     std::vector<gitr_precision> gitrE(nZ,0.0);
    for(int j=0;j<nZ;j++)
    {
      pz[0] = j*dz;
      minDistance =
          getE(px[0], py[0], pz[0], thisE, boundaries.data(), nLines,
               nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
               n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
               &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
               &closeGeom_sheath.front(), closestBoundaryIndex);
      gitrE[j] = thisE[2];
    }
      std::cout << "minDist " << minDistance << std::endl; 
      std::cout << "Efield " << thisE[0] << " " << thisE[1] << " " << thisE[2] << std::endl; 

std::ifstream in( E_FIELD_TEST_FILE );

std::string str;
std::vector<gitr_precision> gold;
// Read the next line from File untill it reaches the end.
while (std::getline(in, str))
{
    // Line contains string of length > 0 then save it in vector
    if(str.size() > 0){
      std::size_t sz = str.length();
      gitr_precision val = std::atof(str.c_str());
        gold.push_back(-val);
        }
}
for(int i=0;i<gold.size();i++){
std::cout << gold[i] << std::endl;
}
gold[0] = 0.0;
    // Compare vectors to ensure reproducibility
    gitr_precision margin = 0.1;
    gitr_precision epsilon = 0.001;
    REQUIRE(compareVectors<gitr_precision>(gitrE,gold,epsilon,margin));
  }
}
