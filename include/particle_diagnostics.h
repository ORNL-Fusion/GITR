#pragma once

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER_DEVICE __device__
#define CUDA_CALLABLE_MEMBER_HOST __host__
#else
#define CUDA_CALLABLE_MEMBER_DEVICE
#define CUDA_CALLABLE_MEMBER_HOST
#endif

#include "Particles.h"
#include "Fields.h"
#include "flags.hpp"
#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/random.h>
#endif

#ifdef __GNUC__
#include <random>
#include <stdlib.h>
#endif

#if USE_DOUBLE
typedef double gitr_precision;
#else
typedef float gitr_precision;
#endif


struct particle_diagnostics {
  Flags *flags;
  Particles *particlesPointer;
  Boundary *boundaryVector;
  bool times_logarithmic;
  gitr_precision bin_edge_0_time;
  gitr_precision bin_edge_1_time;
  gitr_precision bin_edge_dt;
  int n_bins_time;
  gitr_precision *particle_time_histogram;
  bool angle_logarithmic;
  gitr_precision bin_edge_0_angle;
  gitr_precision bin_edge_1_angle;
  gitr_precision bin_edge_dtheta;
  int n_bins_angle;
  gitr_precision *particle_angle_histogram;
  int nSurfaces;
  
  particle_diagnostics(Flags *_flags, Particles *_particlesPointer, Boundary *_boundaryVector, bool _time_logarithmic,
                       gitr_precision _bin_edge_0_time, gitr_precision _bin_edge_1_time, gitr_precision _bin_edge_dt,
                       int _n_bins_time, gitr_precision *_particle_time_histogram, bool _angle_logarithmic, 
                       gitr_precision _bin_edge_0_angle, gitr_precision _bin_edge_1_angle, gitr_precision _bin_edge_dtheta, 
                       int _n_bins_angle, gitr_precision *_particle_angle_histogram, int *_nSurfaces);

  CUDA_CALLABLE_MEMBER_DEVICE
  void operator()(std::size_t indx);
};
