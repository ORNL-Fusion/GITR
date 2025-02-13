#include "boris_data_broker.h"
#include <iostream>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "boris.h"
#include "Surfaces.h"
#include "geometryCheck.h"
#include "boundaryInit.h"
#include "slow_math.h"
#include "constants.h"

std::vector< double > boris_data_broker_0::run_2()
{
  class libconfig_string_query query( GET_E_TEST_FILE );
  class flags f( query );

  int const sheath_efield = 0;
  int const presheath_efield = 1;
  int const surface_potential = 0;
  int const geom_hash_sheath = 0;
  int const use_3d_geom = 0;
  int const cylsymm = 0;

  libconfig::Config cfg_geom;

  cfg_geom.setAutoConvert(true);

  importLibConfig( cfg_geom, GET_E_TEST_FILE );
  int nLines = 1;
  sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

  int nSurfaces = importGeometry( f, cfg_geom, boundaries, cylsymm, 
      surface_potential );

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

  // required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1
  // create a unified setup script
  sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

  // uniform bfield //
  br[ 0 ] = std::cos(M_PI*5.0/180);
  // large bfield in teslas gives smaller gyromotion radius 
  by[ 0 ] = 0;
  bz[ 0 ] = -std::sin(M_PI*5.0/180);;

  // for the uniform efield, set efield to 1000 in z just make the cross product geometry 
  // presheath efield is in the bulk plasma and sheath efield is at the surface of the wall 

  sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
  gitr_precision background_Z = 1;
  gitr_precision background_amu = 2;
  gitr_precision biasPotential = 0;


  std::for_each(boundaries.begin(), boundaries.end() - 1,
      boundary_init(f, background_Z, background_amu, nR_Dens, nZ_Dens,
        DensGridr.data(), DensGridz.data(), ni.data(),
        ne.data(), nR_Bfield, nZ_Bfield,
        bfieldGridr.data(), bfieldGridz.data(), br.data(),
        bz.data(), by.data(), nR_Temp, nZ_Temp,
        TempGridr.data(), TempGridz.data(), ti.data(),
        te.data(), biasPotential, surface_potential,
        cylsymm ));

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
  gitr_precision f_psi = 1.0;
  for(int j=0;j<nZ;j++)
  {
    pz[0] = j*dz;
    minDistance =
      getE(f, px[0], 
           py[0], 
           pz[0], 
           thisE, 
           boundaries.data(), 
           nLines,
           nR_closeGeom_sheath,
           nY_closeGeom_sheath,
           nZ_closeGeom_sheath,
          n_closeGeomElements_sheath,
          &closeGeomGridr_sheath.front(),
          &closeGeomGridy_sheath.front(),
          &closeGeomGridz_sheath.front(),
          &closeGeom_sheath.front(),
          closestBoundaryIndex,
          geom_hash_sheath,
          cylsymm,
          f_psi );

    gitrE[j] = thisE[2];
  }

  std::cout << "minDist " << minDistance << std::endl; 

  std::cout << "Efield " << thisE[0] << " " << thisE[1] << " " << thisE[2] << std::endl; 

  return gitrE;
}

std::vector< double > boris_data_broker_0::run_1()
{
  // Captain! make sure these are consistent with the input file!
  int const sheath_efield = 0;
  int const presheath_efield = 1;
  int const surface_potential = 0;
  int const geom_hash_sheath = 0;
  int const use_3d_geom = 0;
  int const cylsymm = 0;
  // parameters
  int nT = 1e1;

  gitr_precision dt = 1.0e-2;
  gitr_precision vtotal = 1000;

  // total final displacement in y should be vtotal * nT * dt = 1000 * 10 * 1e-2 = 100

  int particle_array_index = 0;

  int initial_x = 0;
  int initial_y = 0;
  int initial_z = 0;

  int velocity_x = 0;
  int velocity_y = vtotal;
  int velocity_z = 0;

  int material_z = 0;

  // amu and error in the "y" dimension appear to be inversely correlated
  gitr_precision amu = 27;

  int charge = 2.0;

  int deprecated_constructor_argument = 0;

  int num_particles = 1;

  // configuration flags
  libconfig::Config cfg_geom;

  cfg_geom.setAutoConvert(true);

  importLibConfig(cfg_geom, BORIS_TEST_FILE_2);

  class libconfig_string_query query_metadata( BORIS_TEST_FILE_2 );
  class flags f( query_metadata );

  // create a particle
  auto particleArray =
    new Particles( num_particles, deprecated_constructor_argument, cfg_geom );

  thrust::counting_iterator<std::size_t> particle_iterator_start(0);

  thrust::counting_iterator<std::size_t> particle_iterator_end(1);

  particleArray->setParticleV( particle_array_index, 
      initial_x,
      initial_y,
      initial_z,
      velocity_x,
      velocity_y,
      velocity_z,
      material_z,
      amu,
      charge,
      dt );

  // dummy variables start

  // hashing dummies
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
  sim::Array<gitr_precision> closeGeomGridr(1);
  sim::Array<gitr_precision> closeGeomGridy(1);
  sim::Array<gitr_precision> closeGeomGridz(1);
  sim::Array<int> closeGeom(1, 0);

  // boundary dummies
  int nLines = 0;
  sim::Array<Boundary> boundaries( nLines + 1, Boundary() );


  int n_closeGeomElements_sheath = 1;

  int nR_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridr_sheath(nR_closeGeom_sheath);

  int nY_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridy_sheath(nY_closeGeom_sheath);

  int nZ_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridz_sheath(nZ_closeGeom_sheath);

  int nGeomHash_sheath = 1;

  sim::Array<int>            closeGeom_sheath(nGeomHash_sheath);

  // presheath efield is in the bulk plasma and sheath efield is at the surface of the wall

  int nR_PreSheathEfield = 1;
  int nY_PreSheathEfield = 1;
  int nZ_PreSheathEfield = 1;
  int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

  int nR_Bfield = 1;
  int nZ_Bfield = 1;
  int n_Bfield = 1;


  // electric field array declarations

  // domain grid
  sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield);
  sim::Array<gitr_precision> preSheathEGridy(nY_PreSheathEfield);
  sim::Array<gitr_precision> preSheathEGridz(nZ_PreSheathEfield);

  // values
  sim::Array<gitr_precision> PSEr(nPSEs); 
  sim::Array<gitr_precision> PSEz(nPSEs); 
  sim::Array<gitr_precision> PSEt(nPSEs);

  // magnetic field array declarations

  // domain grid
  sim::Array<gitr_precision> bfieldGridr(nR_Bfield);
  sim::Array<gitr_precision> bfieldGridz(nZ_Bfield);

  // values
  sim::Array<gitr_precision> br(n_Bfield); 
  sim::Array<gitr_precision> by(n_Bfield);
  sim::Array<gitr_precision> bz(n_Bfield);

  // dummy variables end

  // uniform bfield
  br[ 0 ] = 1;
  by[ 0 ] = 0;
  bz[ 0 ] = 0;

  // r is x 
  // y is t 
  PSEr[ 0 ] = 0;
  PSEz[ 0 ] = 1000;
  PSEt[ 0 ] = 0;


  // create boris operator
  move_boris boris( particleArray,
      dt,
      boundaries.data(),
      nLines,
      nR_Bfield,
      nZ_Bfield,
      bfieldGridr.data(),
      bfieldGridz.data(),
      br.data(),
      bz.data(),
      by.data(),
      nR_PreSheathEfield,
      nY_PreSheathEfield,
      nZ_PreSheathEfield,
      &preSheathEGridr.front(),
      &preSheathEGridy.front(),
      &preSheathEGridz.front(),
      &PSEr.front(),
      &PSEz.front(),
      &PSEt.front(),
      nR_closeGeom_sheath,
      nY_closeGeom_sheath,
      nZ_closeGeom_sheath,
      n_closeGeomElements_sheath,
      closeGeomGridr_sheath.data(),
      &closeGeomGridy_sheath.front(),
      &closeGeomGridz_sheath.front(),
      &closeGeom_sheath.front(),
      f );

  // time loop
  for (int tt = 0; tt < nT; tt++)
  {

    thrust::for_each( thrust::device,
        particle_iterator_start,
        particle_iterator_end,
        boris );


    // manually advance the particle 
    particleArray->xprevious[ 0 ] = particleArray->x[0];
    particleArray->yprevious[ 0 ] = particleArray->y[0];
    particleArray->zprevious[ 0 ] = particleArray->z[0];

  }

  // total final displacement in y should be vtotal * nT * dt = 1000 * 10 * 1e-2 = 100
  std::vector< double > final_position{ particleArray->x[ 0 ],
    particleArray->y[ 0 ],
    particleArray->z[ 0 ] };

  return final_position;
}

void boris_data_broker::run_boris()
{
  thrust::counting_iterator<std::size_t> particle_iterator_start(0);
  thrust::counting_iterator<std::size_t> particle_iterator_end(1);

  for (int i = 0; i < n_timesteps; i++)
  {
    /* save particle velocity/position */
    v_x_test[ i ] = particleArray->vx[ 0 ];
    v_y_test[ i ] = particleArray->vy[ 0 ];
    v_z_test[ i ] = particleArray->vz[ 0 ];

    pos_x_test[ i ] = particleArray->x[ 0 ];
    pos_y_test[ i ] = particleArray->y[ 0 ];
    pos_z_test[ i ] = particleArray->z[ 0 ];

    /* update particle velocity/position */
    thrust::for_each( thrust::device,
        particle_iterator_start,
        particle_iterator_end,
        boris );

    /* manually advance the particle */
    particleArray->xprevious[ 0 ] = particleArray->x[0];
    particleArray->yprevious[ 0 ] = particleArray->y[0];
    particleArray->zprevious[ 0 ] = particleArray->z[0];

  }
}

boris_data_broker::boris_data_broker( Particles *particleArray, 
    int num_particles,
    int n_timesteps,
    double dt,
    class flags &f_init,
    double b_field_x,
    double b_field_y,
    double b_field_z,
    double e_field_x,
    double e_field_y,
    double e_field_z )

  :
n_timesteps( n_timesteps ),
  v_x_test( n_timesteps ),
  v_y_test( n_timesteps ),
  v_z_test( n_timesteps ),
  pos_x_test( n_timesteps ),
  pos_y_test( n_timesteps ),
  pos_z_test( n_timesteps ),
  dt( dt ),
  f( f_init ),
  particleArray( particleArray ),
  num_particles( num_particles ),
  nR_closeGeom(nHashes, 0),
  nY_closeGeom(nHashes, 0),
  nZ_closeGeom(nHashes, 0),
  nHashPoints(nHashes, 0),
  n_closeGeomElements(nHashes, 0),
  closeGeomGridr(1),
  closeGeomGridy(1),
  closeGeomGridz(1),
  closeGeom(1, 0),
  boundaries( nLines + 1, Boundary() ),
  closeGeomGridr_sheath(nR_closeGeom_sheath),
  closeGeomGridy_sheath(nY_closeGeom_sheath),
  closeGeomGridz_sheath(nZ_closeGeom_sheath),
  closeGeom_sheath(nGeomHash_sheath),
  preSheathEGridr(nR_PreSheathEfield),
  preSheathEGridy(nY_PreSheathEfield),
  preSheathEGridz(nZ_PreSheathEfield),
  PSEr(nPSEs),
  PSEz(nPSEs),
  PSEt(nPSEs),
  bfieldGridr(nR_Bfield),
  bfieldGridz(nZ_Bfield),
  br(n_Bfield),
  by(n_Bfield),
  bz(n_Bfield),
  boris( particleArray, // Captain! boris is initialized here!
      dt,
      boundaries.data(),
      nLines,
      nR_Bfield,
      nZ_Bfield,
      bfieldGridr.data(),
      bfieldGridz.data(),
      br.data(),
      bz.data(),
      by.data(),
      nR_PreSheathEfield,
      nY_PreSheathEfield,
      nZ_PreSheathEfield,
      &preSheathEGridr.front(),
      &preSheathEGridy.front(),
      &preSheathEGridz.front(),
      &PSEr.front(),
      &PSEz.front(),
      &PSEt.front(),
      nR_closeGeom_sheath,
      nY_closeGeom_sheath,
      nZ_closeGeom_sheath,
      n_closeGeomElements_sheath,
      closeGeomGridr_sheath.data(),
      &closeGeomGridy_sheath.front(),
      &closeGeomGridz_sheath.front(),
      &closeGeom_sheath.front(),
      f )
{
  /* uniform bfield */
  br[ 0 ] = b_field_x;
  by[ 0 ] = b_field_y;
  bz[ 0 ] = b_field_z;

  /* uniform efield */
  /* r is x */
  /* y is t */
  PSEr[ 0 ] = e_field_x;
  PSEz[ 0 ] = e_field_z;
  PSEt[ 0 ] = e_field_y;

  preSheathEGridr[ 0 ] = 0;
  preSheathEGridr[ 1 ] = -1;
  preSheathEGridy[ 0 ] = 0;
  preSheathEGridy[ 1 ] = -1;
  preSheathEGridz[ 0 ] = 0;
  preSheathEGridz[ 1 ] = -1;

  bfieldGridr[ 0 ] = 0;
  bfieldGridr[ 1 ] = -1;
  bfieldGridz[ 0 ] = 0;
  bfieldGridz[ 1 ] = -1;

}
