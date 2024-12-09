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

#include "boris_data_broker.h"

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
  //SECTION( "compare vx, vy, vz, x, y, z to analytic solution" )
  SECTION( "t0" )
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
      new Particles( num_particles, deprecated_constructor_argument, cfg_geom );

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

    /* time loop */
    boris_data_broker data_broker( particleArray, 
                                   num_particles,
                                   n_timesteps,
                                   dt,
                                   gitr_flags,
                                   sheath_efield,
                                   presheath_efield,
                                   biased_surface,
                                   geom_hash_sheath,
                                   use_3d_geom,
                                   cylsymm,
                                   b_field_x,
                                   b_field_y,
                                   b_field_z,
                                   e_field_x,
                                   e_field_y,
                                   e_field_z );

    data_broker.run_boris();


    /* rmse check the output vs analytical */
    double tolerance = 1e-6;

    REQUIRE( rmse_based_comparison( pos_x, data_broker.pos_x_test, tolerance ) );

    REQUIRE( rmse_based_comparison( pos_y, data_broker.pos_y_test, tolerance ) );

    REQUIRE( rmse_based_comparison( pos_z, data_broker.pos_z_test, tolerance ) );

    REQUIRE( rmse_based_comparison( v_x, data_broker.v_x_test, tolerance ) );

    tolerance = 1e-5;

    REQUIRE( rmse_based_comparison( v_y, data_broker.v_y_test, tolerance ) );

    REQUIRE( rmse_based_comparison( v_z, data_broker.v_z_test, tolerance ) );
  }

  /*

  USE_PRESHEATH_EFIELD=1
  BFIELD_INTERP=1

  */
  //SECTION( "exb drift experiment - Wein Filter" )
  SECTION( "t1" )
  {
    // total final displacement in y should be vtotal * nT * dt = 1000 * 10 * 1e-2 = 100
    
    std::vector< double > gold{ 0, 100, 0 };

    double tolerance = 1e-9;

    boris_data_broker_0 data_broker;

    std::vector< double > final_position = data_broker.run_1();

    REQUIRE( root_mean_squared_error( final_position, gold ) < tolerance );
  }

  //SECTION( "getE tests" )
  SECTION( "t2" )
  {
    int nZ = 10000;

    std::vector<gitr_precision> gitrE(nZ,0.0);

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

    boris_data_broker_0 data_broker;

    std::vector< double > gitr_e = data_broker.run_2();

    REQUIRE(compareVectors<gitr_precision>(gitr_e,gold,epsilon,margin));
  }
}
