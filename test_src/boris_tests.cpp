#include <iostream>
#include <thrust/execution_policy.h>
#include "catch2/catch_all.hpp"
#include "config_interface.h"
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

/*
    REQUIRE(compareVectors<gitr_precision>(gitrE,gold,epsilon,margin));
*/
template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T epsilon, T margin)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) 
  {
    
    bool margin_check = (a[i] != Catch::Approx(b[i]).margin(margin));
    bool epsilon_check = (a[i] != Catch::Approx(b[i]).epsilon(epsilon));
    
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

/*
   test implements prototype script found at:

   "https://github.com/ORNL-Fusion/GITR_prototyping/blob/main/boris_mini_uniform_inputs.py"
*/
/*

  Note: for this test to work, the project must be compiled with 

  GITR_USE_CYLSYMM=0
  GITR_USE_SHEATH_EFIELD=0
  GITR_USE_PRE_SHEATH_EFIELD=1

  This technically will also pass with GITR_USE_PERP_DIFFUSION=2. Unsure which is correct.

 */

TEST_CASE( "Complex Boris Motion" )
{
  int const sheath_efield = 0;
  int const presheath_efield = 1;
  int const biased_surface = 0;
  int const surface_potential = 0;
  int const geom_hash_sheath = 0;
  int const use_3d_geom = 0;
  int const cylsymm = 0;

  /* Testing complex boris motion implemented in the linked script */
  SECTION( "compare vx, vy, vz, x, y, z to analytic solution" )
  {
    /* amu and error in the "y" dimension appear to be inversely correlated */
    int ion_charge = 4.0;

    /* assume tungsten impurities, +4 ion */
    double amu = 184;
    double mass = amu * gitr_constants::dalton;
    double charge = ion_charge * gitr_constants::electron_volt;

    /* generate a magnetic field (b-field) */
    /* units: teslas */
    double b_field_magnitude = 1;

    /* bfield angle 20 degrees in radians */
    double b_field_angle = 20.0 * gitr_constants::pi / 180.0;

    /* instantiate a b-field at the incident angle in the xz plane */
    /* units: teslas */
    double b_field_x = std::cos( b_field_angle ) * b_field_magnitude;
    
    double b_field_y = 0;

    double b_field_z = -1 * std::sin( b_field_angle ) * b_field_magnitude;

    /* rotate coordinates around the y-axis such that z-direction is aligned with the bfield */
    double rotation_angle = gitr_constants::pi / 2 - b_field_angle;

    std::vector< double > const spatial_transform =
    generate_3d_rotation_matrix( rotation_angle );

    std::vector< double > const spatial_inverse_transform =
    generate_3d_rotation_matrix( -1 * rotation_angle );

    std::vector< double > rotated_b_field =
    slow_dot_product( spatial_transform,
                      std::vector< double >{ b_field_x, b_field_y, b_field_z },
                      3,
                      3 );

    /* generate an electric field */
    /* units: V/m */
    double e_field_x = 5e2;
    double e_field_y = 1e3;
    double e_field_z = 1e3;

    std::vector< double > rotated_e_field =
    slow_dot_product( spatial_transform,
                      std::vector< double >{ e_field_x, e_field_y, e_field_z },
                      3,
                      3 );

    /* generate an initial velocity */
    double v_magnitude = 1e4;

    double v0_x = v_magnitude;
    double v0_y = 0;
    double v0_z = 0;

    std::vector< double > rotated_v =
    slow_dot_product( spatial_transform,
                      std::vector< double >{ v0_x, v0_y, v0_z },
                      3,
                      3 );

    /* E cross B drift velocity */
    std::vector< double > v_e_rotated = slow_cross_product( rotated_e_field, rotated_b_field );
    
    std::transform( v_e_rotated.begin(),
                    v_e_rotated.end(),
                    v_e_rotated.begin(), 
                    [ b_field_magnitude ]( double d ) -> double 
                    { return d / ( b_field_magnitude * b_field_magnitude ); } );

    /* next, obtain the component of the rotated velocity orthogonal to the e and b plane: */
    double vsign = rotated_b_field.back() / b_field_magnitude;

    std::vector< double > v_perp( 3 );

    for( int i = 0; i < v_perp.size(); i++ ) v_perp[ i ] = rotated_v[ i ] - v_e_rotated[ i ];

    double norm_v_perp = std::sqrt( slow_dot_product( v_perp, v_perp, 1, v_perp.size() )[ 0 ] );

    /* generate timing */
    double plasma_frequency = b_field_magnitude * charge / mass;

    double const dt = 1 / ( 1e2 * plasma_frequency );

    int n_timesteps = std::ceil( 1e-4 / dt );

    std::cout << "dt: " << dt << " n_timesteps: " << n_timesteps << std::endl;

    /* calculate v_rotated for x, y, and z at each timestep */
    std::vector< double > v_rotated_x( n_timesteps );

    std::vector< double > v_rotated_y( n_timesteps );

    /* z velocity is zero in rotated coordinate system */

    /* calculate position_rotated for x, y, and z at each timestep */
    std::vector< double > pos_rotated_x( n_timesteps );

    std::vector< double > pos_rotated_y( n_timesteps );

    /* z position is zero in rotated coordinate system */

    for( int i = 0; i < n_timesteps; i++ )
    {
      double const t = i * dt;
      double const phase = plasma_frequency * t;

      v_rotated_x[ i ] = -1 * vsign * norm_v_perp * std::cos( phase ) + v_e_rotated[ 0 ];

      v_rotated_y[ i ] = -1 * vsign * norm_v_perp * std::sin( phase ) + v_e_rotated[ 1 ];

      /* z velocity is zero in rotated coordinate system, can be left implicit */

      pos_rotated_x[ i ] =
      norm_v_perp / plasma_frequency * std::sin( phase ) + v_e_rotated[ 0 ] * t;

      pos_rotated_y[ i ] =
      vsign * norm_v_perp / plasma_frequency * ( std::cos( phase ) - 1 ) + v_e_rotated[ 1 ] * t;

      /* z position is zero in rotated coordinate system, can be left implicit */
    }

    /* rotate both position and velocity back into the original space */
    std::vector< double > v_x( n_timesteps );
    std::vector< double > v_y( n_timesteps );
    std::vector< double > v_z( n_timesteps );
    std::vector< double > pos_x( n_timesteps );
    std::vector< double > pos_y( n_timesteps );
    std::vector< double > pos_z( n_timesteps );

    for( int i = 0; i < n_timesteps; i++ )
    {
      std::vector< double > const v_i_xyz =
      slow_dot_product( spatial_inverse_transform,
                        std::vector< double >{ v_rotated_x[ i ], v_rotated_y[ i ], 0 },
                        3,
                        3 );

      std::vector< double > const pos_i_xyz =
      slow_dot_product( spatial_inverse_transform,
                        std::vector< double >{ pos_rotated_x[ i ], pos_rotated_y[ i ], 0 },
                        3,
                        3 );

      v_x[ i ] = v_i_xyz[ 0 ];
      v_y[ i ] = v_i_xyz[ 1 ];
      v_z[ i ] = v_i_xyz[ 2 ];

      pos_x[ i ] = pos_i_xyz[ 0 ];
      pos_y[ i ] = pos_i_xyz[ 1 ];
      pos_z[ i ] = pos_i_xyz[ 2 ];
    }

    /* add electrostatic force on the particle by projecting the e-field onto the b-field and 
       adding the parallel electric force - the perpendicular force was already handled in
       rotated space. This portion is handled in the original space */
    double e_field_parallel = 
    slow_dot_product( std::vector< double >{ e_field_x, e_field_y, e_field_z },
                      std::vector< double >{ b_field_x, b_field_y, b_field_z }, 1, 3 ).front();

    for( int i = 0; i < n_timesteps; i++ )
    {
      double const t = i * dt;

      double const v_parallel = charge * e_field_parallel * t / mass;

      v_x[ i ] += v_parallel * b_field_x / b_field_magnitude;
      v_y[ i ] += v_parallel * b_field_y / b_field_magnitude;
      v_z[ i ] += v_parallel * b_field_z / b_field_magnitude;

      
      pos_x[ i ] += v_parallel * b_field_x / b_field_magnitude * t / 2;
      pos_y[ i ] += v_parallel * b_field_y / b_field_magnitude * t / 2;
      pos_z[ i ] += v_parallel * b_field_z / b_field_magnitude * t / 2;
    }

    /* next, run the GITR boris pusher and get vectors to compare against */
    int particle_array_index = 0;

    int initial_x = 0;
    int initial_y = 0;
    int initial_z = 0;

    /* not sure what the role of this variable is in this context... */
    int material_z = 0;

    int deprecated_constructor_argument = 0;

    int num_particles = 1;

    /* configuration flags */
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, BORIS_TEST_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    /* create a particle */
    auto particleArray =
      new Particles( num_particles, deprecated_constructor_argument, cfg_geom, gitr_flags );

    thrust::counting_iterator<std::size_t> particle_iterator_start(0);

    thrust::counting_iterator<std::size_t> particle_iterator_end(1);

    particleArray->setParticleV( particle_array_index, 
        initial_x,
        initial_y,
        initial_z,
        v_x[ 0 ],
        v_y[ 0 ],
        v_z[ 0 ],
        material_z,
        amu,
        ion_charge,
        dt );

    /* dummy variables start */

    /* hashing dummies */
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

    /* boundary dummies */
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

    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    int nR_PreSheathEfield = 1;
    int nY_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;
    int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    int n_Bfield = 1;


    /* electric field array declarations */

    /* domain grid */
    sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield);
    sim::Array<gitr_precision> preSheathEGridy(nY_PreSheathEfield);
    sim::Array<gitr_precision> preSheathEGridz(nZ_PreSheathEfield);

    /* values */
    sim::Array<gitr_precision> PSEr(nPSEs); 
    sim::Array<gitr_precision> PSEz(nPSEs); 
    sim::Array<gitr_precision> PSEt(nPSEs);

    /* magnetic field array declarations */

    /* domain grid */
    sim::Array<gitr_precision> bfieldGridr(nR_Bfield);
    sim::Array<gitr_precision> bfieldGridz(nZ_Bfield);

    /* values */
    sim::Array<gitr_precision> br(n_Bfield); 
    sim::Array<gitr_precision> by(n_Bfield);
    sim::Array<gitr_precision> bz(n_Bfield);

    /* dummy variables end */

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


    /* create boris operator */
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
        gitr_flags,
        sheath_efield,
        presheath_efield,
        biased_surface,
        geom_hash_sheath,
        use_3d_geom,
        cylsymm );

    /* time loop */
    std::vector< double > v_x_test( n_timesteps );
    std::vector< double > v_y_test( n_timesteps );
    std::vector< double > v_z_test( n_timesteps );

    std::vector< double > pos_x_test( n_timesteps );
    std::vector< double > pos_y_test( n_timesteps );
    std::vector< double > pos_z_test( n_timesteps );

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

    /* rmse check the output vs analytical */
    double tolerance = 1e-6;

    REQUIRE( rmse_based_comparison( pos_x, pos_x_test, tolerance ) );

    REQUIRE( rmse_based_comparison( pos_y, pos_y_test, tolerance ) );

    REQUIRE( rmse_based_comparison( pos_z, pos_z_test, tolerance ) );

    REQUIRE( rmse_based_comparison( v_x, v_x_test, tolerance ) );

    tolerance = 1e-5;

    REQUIRE( rmse_based_comparison( v_y, v_y_test, tolerance ) );

    REQUIRE( rmse_based_comparison( v_z, v_z_test, tolerance ) );
  }

  /*

  USE_PRESHEATH_EFIELD=1
  BFIELD_INTERP=1

  */
  SECTION( "exb drift experiment - Wein Filter" )
  {
    /* parameters */
    int nT = 1e1;

    gitr_precision dt = 1.0e-2;
    gitr_precision vtotal = 1000;

    /* total final displacement in y should be vtotal * nT * dt = 1000 * 10 * 1e-2 = 100 */

    int particle_array_index = 0;

    int initial_x = 0;
    int initial_y = 0;
    int initial_z = 0;

    int velocity_x = 0;
    int velocity_y = vtotal;
    int velocity_z = 0;

    int material_z = 0;

    /* amu and error in the "y" dimension appear to be inversely correlated */
    gitr_precision amu = 27;

    int charge = 2.0;

    int deprecated_constructor_argument = 0;

    int num_particles = 1;

    /* configuration flags */
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, BORIS_TEST_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    /* create a particle */
    auto particleArray =
    new Particles( num_particles, deprecated_constructor_argument, cfg_geom, gitr_flags );

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

    /* dummy variables start */

    /* hashing dummies */
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

    /* boundary dummies */
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

    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    int nR_PreSheathEfield = 1;
    int nY_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;
    int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

    int nR_Bfield = 1;
    int nZ_Bfield = 1;
    int n_Bfield = 1;


    /* electric field array declarations */

    /* domain grid */
    sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield);
    sim::Array<gitr_precision> preSheathEGridy(nY_PreSheathEfield);
    sim::Array<gitr_precision> preSheathEGridz(nZ_PreSheathEfield);

    /* values */
    sim::Array<gitr_precision> PSEr(nPSEs); 
    sim::Array<gitr_precision> PSEz(nPSEs); 
    sim::Array<gitr_precision> PSEt(nPSEs);

    /* magnetic field array declarations */
    
    /* domain grid */
    sim::Array<gitr_precision> bfieldGridr(nR_Bfield);
    sim::Array<gitr_precision> bfieldGridz(nZ_Bfield);

    /* values */
    sim::Array<gitr_precision> br(n_Bfield); 
    sim::Array<gitr_precision> by(n_Bfield);
    sim::Array<gitr_precision> bz(n_Bfield);

    /* dummy variables end */

    /* uniform bfield */
    br[ 0 ] = 1;
    by[ 0 ] = 0;
    bz[ 0 ] = 0;

    /* r is x */
    /* y is t */
    PSEr[ 0 ] = 0;
    PSEz[ 0 ] = 1000;
    PSEt[ 0 ] = 0;


    /* create boris operator */
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
                      gitr_flags,
                      sheath_efield,
                      presheath_efield,
                      biased_surface,
                      geom_hash_sheath,
                      use_3d_geom,
                      cylsymm );

    /* time loop */
    for (int tt = 0; tt < nT; tt++)
    {

      thrust::for_each( thrust::device,
                        particle_iterator_start,
                        particle_iterator_end,
                        boris );


      /* manually advance the particle */
      particleArray->xprevious[ 0 ] = particleArray->x[0];
      particleArray->yprevious[ 0 ] = particleArray->y[0];
      particleArray->zprevious[ 0 ] = particleArray->z[0];

    }

    /* total final displacement in y should be vtotal * nT * dt = 1000 * 10 * 1e-2 = 100 */
    std::vector< double > final_position{ particleArray->x[ 0 ],
                                          particleArray->y[ 0 ],
                                          particleArray->z[ 0 ] };
    
    std::vector< double > gold{ 0, 100, 0 };

    double tolerance = 1e-9;

    REQUIRE( root_mean_squared_error( final_position, gold ) < tolerance );
  }

  SECTION( "getE tests" )
  {
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig( cfg_geom, GET_E_TEST_FILE );
    int nLines = 1;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries, use_3d_geom, cylsymm, 
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
                              te.data(), biasPotential, biased_surface, surface_potential,
                              use_3d_geom, cylsymm ));
    
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
               &closeGeom_sheath.front(), closestBoundaryIndex, biased_surface,
               use_3d_geom, geom_hash_sheath, cylsymm );
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
        if(str.size() > 0)
        {
          gitr_precision val = std::atof(str.c_str());
          
          gold.push_back(-val);
        }
    }

    for(int i=0;i<gold.size();i++)
    {
      std::cout << gold[i] << " " << gitrE[ i ] << std::endl;
    }

    gold[0] = 0.0;

    // Compare vectors to ensure reproducibility
    gitr_precision margin = 0.1;
    gitr_precision epsilon = 0.001;
    REQUIRE(compareVectors<gitr_precision>(gitrE,gold,epsilon,margin));
  }
}




















