#include "array.h"
#include "boris.h"
#include "Particles.h"
#include "flags.hpp"
#include <vector>
#include "config_interface.h"

class boris_data_broker_0
{
  public:

  boris_data_broker_0(){}

  std::vector< double > run_2();

  std::vector< double > run_1();

  std::vector< double > run_0();
};

class boris_data_broker
{
  public:

  /* default constructor */
  boris_data_broker( Particles *particleArray,
                     int num_particles,
                     int n_timesteps,
                     double dt,
                     Flags *flags,
                     class flags &f_init,
                     int sheath_efield,
                     int presheath_efield,
                     int biased_surface,
                     int geom_hash_sheath,
                     int use_3d_geom,
                     int cylsymm,
                     double b_field_x,
                     double b_field_y,
                     double b_field_z,
                     double e_field_x,
                     double e_field_y,
                     double e_field_z );

  /* Carries out the transformation */
  void run_boris();

  /* Wait a minute... How can these even be here?? */
  std::vector< double > v_x_test;//( n_timesteps );
  std::vector< double > v_y_test;//( n_timesteps );
  std::vector< double > v_z_test;//( n_timesteps );

  std::vector< double > pos_x_test;//( n_timesteps );
  std::vector< double > pos_y_test;//( n_timesteps );
  std::vector< double > pos_z_test;//( n_timesteps );

  Particles *particleArray;

  int num_particles;

  int n_timesteps;

  double dt;

  Flags *flags;
  class flags &f;

  /* Captain! Are these allocated on device? */
  /* hashing dummies */
  int nHashes = 1;
  sim::Array<int> nR_closeGeom;//(nHashes, 0);
  sim::Array<int> nY_closeGeom;//(nHashes, 0);
  sim::Array<int> nZ_closeGeom;//(nHashes, 0);
  sim::Array<int> nHashPoints;//(nHashes, 0);
  sim::Array<int> n_closeGeomElements;//(nHashes, 0);
  int nEdist = 1;
  gitr_precision E0dist = 0.0;
  gitr_precision Edist = 0.0;
  int nAdist = 1;
  gitr_precision A0dist = 0.0;
  gitr_precision Adist = 0.0;
  sim::Array<gitr_precision> closeGeomGridr;//(1);
  sim::Array<gitr_precision> closeGeomGridy;//(1);
  sim::Array<gitr_precision> closeGeomGridz;//(1);
  sim::Array<int> closeGeom;//(1, 0);

  /* boundary dummies */
  int nLines = 0;
  sim::Array<Boundary> boundaries;//( nLines + 1, Boundary() );

  int n_closeGeomElements_sheath = 1;

  int nR_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridr_sheath;//(nR_closeGeom_sheath);

  int nY_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridy_sheath;//(nY_closeGeom_sheath);

  int nZ_closeGeom_sheath = 1;

  sim::Array<gitr_precision> closeGeomGridz_sheath;//(nZ_closeGeom_sheath);

  int nGeomHash_sheath = 1;

  sim::Array<int>            closeGeom_sheath;//(nGeomHash_sheath);

  /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

  int const nR_PreSheathEfield = 1;
  int const nY_PreSheathEfield = 1;
  int const nZ_PreSheathEfield = 1;
  int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

  int const nR_Bfield = 1;
  int const nZ_Bfield = 1;
  int n_Bfield = 1;


  /* electric field array declarations */

  /* domain grid */
  sim::Array<gitr_precision> preSheathEGridr;//(nR_PreSheathEfield);
  sim::Array<gitr_precision> preSheathEGridy;//(nY_PreSheathEfield);
  sim::Array<gitr_precision> preSheathEGridz;//(nZ_PreSheathEfield);

  /* values */
  sim::Array<gitr_precision> PSEr;//(nPSEs); 
  sim::Array<gitr_precision> PSEz;//(nPSEs); 
  sim::Array<gitr_precision> PSEt;//(nPSEs);

  /* magnetic field array declarations */
  
  /* domain grid */
  sim::Array<gitr_precision> bfieldGridr;//(nR_Bfield);
  sim::Array<gitr_precision> bfieldGridz;//(nZ_Bfield);

  /* values */
  sim::Array<gitr_precision> br;//(n_Bfield); 
  sim::Array<gitr_precision> by;//(n_Bfield);
  sim::Array<gitr_precision> bz;//(n_Bfield);

  move_boris boris;
};
