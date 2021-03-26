#ifndef _FIELD_
#define _FIELD_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#if USE_MPI > 0
#include "mpi.h"
#include "netcdf.h"
#include "netcdf_par.h"
#endif
#include "utils.h"
#include <cstdlib>
#include <libconfig.h++>
#include <netcdf>
#include <stdio.h>
#include <vector>
#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#endif

#include <random>
enum FieldType { FieldType_constant = 0, FieldType_2d_xz = 1, FieldType_2d_rz = 2 };
template <typename T> using FunctionHandler = float (T::*)(float, float, float);

class Field : public ManagedAllocation {

public:
  FieldType field_type;
  std::size_t nD;
  sim::Array<std::size_t> dimensions;
  sim::Array<float> x;
  sim::Array<float> y;
  sim::Array<float> z;
  sim::Array<float> values;
  FunctionHandler<Field> function;

  CUDA_CALLABLE_MEMBER Field();

  CUDA_CALLABLE_MEMBER Field(libconfig::Config &cfg, std::string field_name); 
  
  CUDA_CALLABLE_MEMBER float return_const(float x, float y, float z);

  CUDA_CALLABLE_MEMBER float returnOne(float x, float y, float z);
  
  CUDA_CALLABLE_MEMBER float cylindrical_2d_interpolation(float x, float y, float z);

  float interp_2d(float x, float y); 

  CUDA_CALLABLE_MEMBER FunctionHandler<Field> returnPointerTable(int num);

  float interpolate(float x, float y, float z); 
};

#endif
