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
#include <netcdf_par.h>
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

  CUDA_CALLABLE_MEMBER
  Field()
      : field_type{FieldType_2d_xz}, nD{0}, dimensions{3, 0}, x{2, 0.0}, y{3,
                                                                           0.0},
        z{3, 0.0}, values{9, 0.0}, function{returnPointerTable(1)} {};

  CUDA_CALLABLE_MEMBER
  Field(libconfig::Config &cfg, std::string field_name) {
    int interpolator_number = get_variable<int>(cfg, field_name + ".interpolation");
    if ( interpolator_number == 0)
    {
      //std::cout << "constant " << std::endl;
      field_type = FieldType_constant;
      nD = 0;
      values[0] = get_variable<float>(cfg, field_name+ ".value");
      function = returnPointerTable(field_type);
    }
    else if (interpolator_number > 0)
    {
      //std::cout << "non - constant " << std::endl;
    std::string ncfilename =
        get_variable<const char *>(cfg, field_name + ".filename");
    std::string ncvarname =
        get_variable<const char *>(cfg, field_name + ".variable_name");

    int err, status, ncid, groupid, cmode, val_id, r_id, retval;
    cmode = NC_NOWRITE;
#if USE_MPI == 1
    err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL,
                      &ncid);
#else
    err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
    err = nc_inq_varid(ncid, ncvarname.c_str(), &val_id);

    // get number of dimensions
    int ndimsp;
    nc_inq_varndims(ncid, val_id, &ndimsp);
    nD = ndimsp;

    // get dimension lengths
    dimensions.resize(nD);
    int *dimidsp = new int[ndimsp];
    nc_inq_vardimid(ncid, val_id, &dimidsp[0]);

    std::vector<std::string> dimension_names(ndimsp);
    std::vector<std::size_t> dimension_lengths(ndimsp);

    // get dimension names and lengths and compute
    // value length
    int val_size = 1;
    for (int i = 0; i < ndimsp; i++) {
      char *name = new char[NC_MAX_NAME + 1];
      err = nc_inq_dimname(ncid, i, name);
      std::string name_string(name);
      dimension_names[i] = name_string;
      err = nc_inq_dimlen(ncid, i, &dimension_lengths[i]);
      val_size = val_size * dimension_lengths[i];
    }
    dimensions = dimension_lengths;

    // Determine field type from dimension names
    const std::string char_r = "r";
    const std::string char_z = "z";
    if ((nD == 2) & (dimension_names[0] == char_r) & (dimension_names[1] == char_z)) {
      field_type = FieldType_2d_rz;
    }
    function = returnPointerTable(field_type);

    // Fill grids and values of vectors
    for (int index = 0; index < nD; index++) {
      if (index <= dimension_lengths.size() - 1 || index == - 1) {
        int vec_size;
        std::string var_name;
        var_name = dimension_names[index];
        vec_size = dimension_lengths[index];
        int val_id;
        err = nc_inq_varid(ncid, var_name.c_str(), &val_id);
        std::vector<float> vec(vec_size, 0.0);
        err = nc_get_var_float(ncid, val_id, &vec[0]);

        if (index == 0) {
          x.resize(vec_size);
          x = vec;
        } 
	else if (index == 1) {
          y.resize(vec_size);
          y = vec;
        } 
	else if (index == 2) {
          z.resize(vec_size);
          z = vec;
        }

      } else {
        std::vector<float> vec0(1, 0.0);
      }
    }

    // get values
    err = nc_inq_varid(ncid, ncvarname.c_str(), &val_id);
    std::vector<float> vec(val_size, 0.0);
    err = nc_get_var_float(ncid, val_id, &vec[0]);
    values.resize(val_size);
    values = vec;
    nc_close(ncid);
   }
  };
  
  
  CUDA_CALLABLE_MEMBER
  float return_const(float x, float y, float z) { return values[0]; }

  CUDA_CALLABLE_MEMBER
  float returnOne(float x, float y, float z) { return 1; }
  
  CUDA_CALLABLE_MEMBER
  float cylindrical_2d_interpolation(float x, float y, float z) {
    float r = sqrt(x * x + y * y);
    float value = interp_2d(r, z);
    return value;
  }

  // CUDA_CALLABLE_MEMBER
  float interp_2d(float x, float y) {
    float dx = this->x[1] - this->x[0];
    float dy = this->y[1] - this->y[0];
    
    int nx = this->dimensions[0];
    int ny = this->dimensions[1];

    int x_index = std::floor((x - this->x[0]) / dx);
    int y_index = std::floor((y - this->y[0]) / dy);
    
    float value_y0 =
        ((this->x[x_index + 1] - x) * this->values[x_index * ny + y_index] +
         (x - this->x[x_index]) * this->values[(x_index + 1) * ny + y_index]) /
        dx;

    float value_y1 = ((this->x[x_index + 1] - x) *
                          this->values[x_index * ny + (y_index + 1)] +
                      (x - this->x[x_index]) *
                          this->values[(x_index + 1) * ny + (y_index + 1)]) /
                     dx;
    
    float value = ((this->y[y_index + 1] - y) * value_y0 +
                   (y - this->y[y_index]) * value_y1) /
                  dy;
    return value;
  }

  CUDA_CALLABLE_MEMBER
  FunctionHandler<Field> returnPointerTable(int num) {
    if (num == FieldType_constant) {
      return &Field::return_const;
    }
    else if(num == FieldType_2d_xz) {
      return &Field::returnOne;
    } else {
      return &Field::cylindrical_2d_interpolation;
    }
  }

  float interpolate(float x, float y, float z) {
    float val = (this->*(this->function))(x, y, z);
    return val;
  }

};

#endif
