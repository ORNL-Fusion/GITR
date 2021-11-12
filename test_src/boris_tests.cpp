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
#include "slow_math.h"

/* today: boris test, wein filter test, fix catch stuff, slow math implementations */

/* test implements prototype script: */
/*
   "https://github.com/ORNL-Fusion/GITR_prototyping/blob/main/boris_mini_uniform_inputs.py"
*/
TEST_CASE( "Complex Boris Motion" )
{
  /* Testing complex boris motion implemented in the linked script */
  SECTION( "t0" )
  {
    //double const pi = 3.14159265358979323846;
    double mass_proton = 1.6605e-27;
    double mass_electron = 9.109383701528e-31;
    double charge_electron = 1.602176634e-19;

    /* assume tungsten impurities, +4 ion */
    double mass = 184 * mass_proton;
    double charge = 4 * charge_electron;

    /* generate a magnetic field (b-field) */
    /* units: teslas */
    double b_field_magnitude = 0.5;

    /* bfield angle 20 degrees in radians */
    double b_field_angle = 20.0 * pi / 180.0;

    /* instantiate a b-field at the incident angle in the xz plane */
    /* units: teslas */
    double b_field_x = std::cos( b_field_angle ) * b_field_magnitude;
    
    double b_field_y = 0;

    double b_field_z = -1 * std::sin( b_field_angle ) * b_field_magnitude;

    /* rotate coordinates around the y-axis such that z-direction is aligned with the bfield */
    double rotation_angle = pi / 2 - b_field_angle;

    std::vector< double > const spatial_transform =
    generate_3d_rotation_matrix( rotation_angle );

    std::vector< double > const spatial_inverse_transform =
    generate_3d_rotation_matrix( -1 * rotation_angle );

    std::vector< double > rotated_b_field =
    slow_dot_product( spatial_transform,
                      std::vector< double >{ b_field_x, b_field_y, b_field_z },
                      3,
                      3 );

    /* check here - seems correct */
    std::cout << "Ahoy, Captain! bfield: "
              << b_field_x << " " << b_field_y << " " << b_field_z << std::endl;

    std::cout << "Ahoy, Captain! transformed bfield: " 
              << rotated_b_field[ 0 ]
              << " "
              << rotated_b_field[ 1 ]
              << " "
              << rotated_b_field[ 2 ]
              << std::endl;

    /* Captain! Here */
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

    std::cout << "Ahoy, Captain! transformed efield: " 
              << rotated_e_field[ 0 ]
              << " "
              << rotated_e_field[ 1 ]
              << " "
              << rotated_e_field[ 2 ]
              << std::endl;

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

    std::cout << "Ahoy, Captain! transformed v_e_rotated: " 
              << v_e_rotated[ 0 ]
              << " "
              << v_e_rotated[ 1 ]
              << " "
              << v_e_rotated[ 2 ]
              << std::endl;

    /* next, obtain the component of the rotated velocity orthogonal to the e and b plane: */
    double vsign = rotated_b_field.back() / b_field_magnitude;

    std::vector< double > v_perp( 3 );

    for( int i = 0; i < v_perp.size(); i++ ) v_perp[ i ] = rotated_v[ i ] - v_e_rotated[ i ];

    double norm_v_perp = std::sqrt( slow_dot_product( v_perp, v_perp, 1, v_perp.size() )[ 0 ] );

    std::cout << "Ahoy, Captain! norm_v_perp: " << norm_v_perp << std::endl;

    /* generate timing */
    double plasma_frequency = b_field_magnitude * charge / mass;

    double const dt = 1 / ( 1e4 * plasma_frequency );
    int n_timesteps = std::ceil( 1e-4 / dt );

    std::cout << "Ahoy, Captain! plasma_frequency: " << plasma_frequency << std::endl;
    std::cout << "Ahoy, Captain! dt: " << dt << std::endl;
    std::cout << "Ahoy, Captain! n_timesteps: " << n_timesteps << std::endl;

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

      /* Captain! Question for tomorrow: why is vsign not on the first one here? */
      pos_rotated_x[ i ] =
      norm_v_perp / plasma_frequency * std::sin( phase ) + v_e_rotated[ 0 ] * t;

      pos_rotated_y[ i ] =
      vsign * norm_v_perp / plasma_frequency * ( std::cos( phase ) - 1 ) + v_e_rotated[ 1 ] * t;

      /* z position is zero in rotated coordinate system, can be left implicit */
    }

    /* spot check */
    std::cout << "Ahoy, Captain! v_rotated_x[ 500000 ] " << v_rotated_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! v_rotated_y[ 500000 ] " << v_rotated_y[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_rotated_x[ 500000 ] " << pos_rotated_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_rotated_y[ 500000 ] " << pos_rotated_y[ 500000 ] << std::endl;
    
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

    /* spot check */
    std::cout << "Ahoy, Captain! v_x[ 500000 ] " << v_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! v_y[ 500000 ] " << v_y[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! v_z[ 500000 ] " << v_z[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_x[ 500000 ] " << pos_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_y[ 500000 ] " << pos_y[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_z[ 500000 ] " << pos_z[ 500000 ] << std::endl;

    /* add electrostatic force on the particle by project the e-field onto the b-field and 
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

    /* spot check again */
    std::cout << "After factoring in parallel electrostatic force:" << std::endl;
    std::cout << "Ahoy, Captain! v_x[ 500000 ] " << v_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! v_y[ 500000 ] " << v_y[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! v_z[ 500000 ] " << v_z[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_x[ 500000 ] " << pos_x[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_y[ 500000 ] " << pos_y[ 500000 ] << std::endl;
    std::cout << "Ahoy, Captain! pos_z[ 500000 ] " << pos_z[ 500000 ] << std::endl;

    /* finished with analytical solution! */

    /* next, run the GITR boris pusher and get vectors to compare against */


    /* rmse check the output vs analytical */

    REQUIRE( false );
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
                      gitr_flags );

    /* get particle xyz before */
    /* time loop */
    std::cout << "Before: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;

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

    std::cout << "After: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;
    /* get particle xyz after */

    /* what should it be analytically? */

    /* charge coulombs, e in V/m, b in Teslas */
    /* q * e / b^2 */
  }
}
