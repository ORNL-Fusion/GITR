#ifndef _BORIS_
#define _BORIS_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#include "thrust/extrema.h"
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_CALLABLE_MEMBER_DEVICE
#endif

#include <algorithm>
#include "Particles.h"
#include "Boundary.h"
#include "interp2d.hpp"
#include "utils.h"
#include <algorithm>
#include "interpolator.h"
#include "config_interface.h"

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif

CUDA_CALLABLE_MEMBER
void vectorAdd(gitr_precision A[], gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorSubtract(gitr_precision A[], gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorScalarMult(gitr_precision a, gitr_precision B[],gitr_precision C[]);

CUDA_CALLABLE_MEMBER
void vectorAssign(gitr_precision a, gitr_precision b,gitr_precision c, gitr_precision D[]);

CUDA_CALLABLE_MEMBER
gitr_precision vectorNorm(gitr_precision A[]);

CUDA_CALLABLE_MEMBER
void vectorNormalize(gitr_precision A[],gitr_precision B[]);

CUDA_CALLABLE_MEMBER
gitr_precision vectorDotProduct(gitr_precision A[], gitr_precision B[]);

CUDA_CALLABLE_MEMBER
void vectorCrossProduct(gitr_precision A[], gitr_precision B[], gitr_precision C[]);

CUDA_CALLABLE_MEMBER
gitr_precision getE ( class flags config_flags, gitr_precision x0, gitr_precision y, gitr_precision z, gitr_precision E[], Boundary *boundaryVector, int nLines,
       int nR_closeGeom, int nY_closeGeom,int nZ_closeGeom, int n_closeGeomElements, 
       gitr_precision *closeGeomGridr,gitr_precision *closeGeomGridy, gitr_precision *closeGeomGridz, int *closeGeom, 
         int&  closestBoundaryIndex,
         gitr_precision& f_psi ); 

struct move_boris { 

    /* Ahoy, Captain! new code for class based boris interpolation */
    CUDA_CALLABLE_MEMBER    
    gitr_precision interp2dCombined ( gitr_precision x,
                                      gitr_precision y,
                                      gitr_precision z,
                                      int nx,
                                      int nz,
                                      gitr_precision* gridx,
                                      gitr_precision* gridz,
                                      gitr_precision* data,
                                      int cylsymm );

    CUDA_CALLABLE_MEMBER    
    void interp2dVector( gitr_precision* field, 
                         gitr_precision x,
                         gitr_precision y,
                         gitr_precision z,
                         int nx,
                         int nz,
                         gitr_precision* gridx,
                         gitr_precision* gridz,
                         gitr_precision* datar,
                         gitr_precision* dataz,
                         gitr_precision* datat,
                         int cylsymm );

    Particles *particlesPointer;
    //int& tt;
    Boundary *boundaryVector;
    int nR_Bfield;
    int nZ_Bfield;
    gitr_precision * BfieldGridRDevicePointer;
    gitr_precision * BfieldGridZDevicePointer;
    gitr_precision * BfieldRDevicePointer;
    gitr_precision * BfieldZDevicePointer;
    gitr_precision * BfieldTDevicePointer;
    int nR_Efield;
    int nY_Efield;
    int nZ_Efield;
    gitr_precision * EfieldGridRDevicePointer;
    gitr_precision * EfieldGridYDevicePointer;
    gitr_precision * EfieldGridZDevicePointer;
    gitr_precision * EfieldRDevicePointer;
    gitr_precision * EfieldZDevicePointer;
    gitr_precision * EfieldTDevicePointer;
    int nR_closeGeom_sheath;
    int nY_closeGeom_sheath;
    int nZ_closeGeom_sheath;
    int n_closeGeomElements_sheath;
    gitr_precision* closeGeomGridr_sheath;
    gitr_precision* closeGeomGridy_sheath;
    gitr_precision* closeGeomGridz_sheath;
    int* closeGeom_sheath; 
    flags config_flags;
    gitr_precision max_dt;

    const gitr_precision span;
    const int nLines;
    gitr_precision magneticForce[3];
    gitr_precision electricForce[3];
    int sheath_efield;
    int presheath_efield;
    int geom_hash_sheath;
    dummy< gitr_precision > dum;
    /* what is a proof of concept experiment you could do? */
    /* what is the next step after that? */
    /*
        double dims[ 2 ];
        double lower_bounds[2];
        double upper_bounds[2];
    interpolated_field< gitr_precision > 
    ifield;
    */

    move_boris(Particles *_particlesPointer, gitr_precision _span, Boundary *_boundaryVector,int _nLines,
            int _nR_Bfield, int _nZ_Bfield,
            gitr_precision * _BfieldGridRDevicePointer,
            gitr_precision * _BfieldGridZDevicePointer,
            gitr_precision * _BfieldRDevicePointer,
            gitr_precision * _BfieldZDevicePointer,
            gitr_precision * _BfieldTDevicePointer,
            int _nR_Efield,int _nY_Efield, int _nZ_Efield,
            gitr_precision * _EfieldGridRDevicePointer,
            gitr_precision * _EfieldGridYDevicePointer,
            gitr_precision * _EfieldGridZDevicePointer,
            gitr_precision * _EfieldRDevicePointer,
            gitr_precision * _EfieldZDevicePointer,
            gitr_precision * _EfieldTDevicePointer,
            int _nR_closeGeom, int _nY_closeGeom,int _nZ_closeGeom, 
            int _n_closeGeomElements, gitr_precision *_closeGeomGridr,
            gitr_precision *_closeGeomGridy, gitr_precision *_closeGeomGridz, 
            int *_closeGeom, 
            flags &f_init,
            gitr_precision _max_dt = 1.0e5);

    CUDA_CALLABLE_MEMBER    
    void operator()(std::size_t indx); 
};

#endif
