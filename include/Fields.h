#ifndef _FIELD_
#define _FIELD_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <libconfig.h++>
#include "mpi.h"
#include "utils.h"
#include <netcdf>
#include <netcdf_par.h>
#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#endif

#include <random>
enum FieldType { 
    FieldType_2d_xz = 1,    FieldType_2d_rz  = 2
};
template<typename T> using FunctionHandler = float(T::*)(float,float,float);

class Field : public ManagedAllocation {
public:
  FieldType field_type;
  std::size_t nD;
  sim::Array<std::size_t> dimensions;
  sim::Array<float> x;
  sim::Array<float> y;
  sim::Array<float> z;
  sim::Array<float> values;
  FunctionHandler<Field> fooHandler;
 //static Field* Create(FieldType type);
 
  CUDA_CALLABLE_MEMBER 
  Field() :
    field_type{FieldType_2d_xz},	  
    nD{0}, dimensions{3,0},
    x{2,0.0},y{3,0.0},z{3,0.0},values{9,0.0},
   fooHandler{returnPointerTable(1)} 
    {};
  
  //CUDA_CALLABLE_MEMBER 
  //Field ( libconfig::Config &cfg, std::string field_name) :
  //  field_type{get_field_type(cfg, field_name)},	  
  //  nD{get_field_n_dimensions(cfg, field_name)},
  //  dimensions{get_field_dimensions(cfg, field_name)},
  //  dimension_names{get_field_dimension_names(cfg, field_name)},
  //  x{get_field_vector(cfg, field_name, 0)},
  //  y{get_field_vector(cfg, field_name, 1)},
  //  z{get_field_vector(cfg, field_name, 2)},
  //  values{get_field_vector(cfg, field_name, -1)},
  // fooHandler{returnPointerTable(field_type)} 
  //  {};
  
  CUDA_CALLABLE_MEMBER 
  Field ( libconfig::Config &cfg, std::string field_name) //:
    //field_type{get_field_type(cfg, field_name)},	  
    //nD{get_field_n_dimensions(cfg, field_name)},
    //dimensions{get_field_dimensions(cfg, field_name)},
    //dimension_names{get_field_dimension_names(cfg, field_name)},
    //x{get_field_vector(cfg, field_name, 0)},
    //y{get_field_vector(cfg, field_name, 1)},
    //z{get_field_vector(cfg, field_name, 2)},
    //values{get_field_vector(cfg, field_name, -1)},
    //fooHandler{returnPointerTable(field_type)} 
    {
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      
      // get number of dimensions
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      nD = ndimsp;

      // get dimension lengths
      dimensions.resize(nD);
      int* dimidsp = new int[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::vector<std::string> dimension_names(ndimsp); 
      std::vector<std::size_t> dimension_lengths(ndimsp); 
      int val_size = 1;
      for(int i=0; i<ndimsp; i++)
      {
	      char* name = new char[NC_MAX_NAME+1];
              err = nc_inq_dimname(ncid,i, name);
	      std::string name_string(name);
	      dimension_names[i] = name_string;
              err = nc_inq_dimlen(ncid,i, &dimension_lengths[i]);
	val_size = val_size*dimension_lengths[i];
      }
      dimensions = dimension_lengths;
     
     // Determine field type from dimension names
      const std::string char_r = "r";
      const std::string char_z = "z";
      if(nD == 2 & dimension_names[0] == char_r
                  & dimension_names[1] == char_z)
      {
       field_type = FieldType_2d_rz;
      }
     fooHandler = returnPointerTable(field_type); 

     // Fill grids and values of vectors
     for(int index=0; index<nD; index++){
     if (index <= dimension_lengths.size() - 1 || index == -1 )
     {
       int vec_size;     
       std::string var_name;
	       var_name = dimension_names[index];   
	       vec_size = dimension_lengths[index];
       int val_id;
       err=nc_inq_varid (ncid, var_name.c_str(), &val_id);
       std::vector<float> vec(vec_size,0.0);
       err = nc_get_var_float( ncid, val_id, &vec[0]);
       if(index ==0){
	       x.resize(vec_size);
	       x = vec;
       }
       else if(index ==1){
	       y.resize(vec_size);
	       y = vec;
       }
       else if(index ==2){
	       z.resize(vec_size);
	       z = vec;
       }

     }
     else{
	     std::vector<float> vec0(1,0.0);
     }
    } 
       
	       //std::cout << "inside -1 " << std::endl;
       err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
       std::vector<float> vec(val_size,0.0);
       err = nc_get_var_float( ncid, val_id, &vec[0]);
       values.resize(val_size);
       values = vec;
       nc_close(ncid);
	  //nD = get_field_n_dimensions(cfg, field_name);
	  //dimensions.resize(nD);
	  //dimensions = get_field_dimensions(cfg, field_name);
    };
//  CUDA_CALLABLE_MEMBER 
//  Field ( libconfig::Config &cfg, std::string field_name) : dimensions{0,0}, x{0,0.0}, y{0,0.0}, z{0,0.0}, values{0,0.0} {
//      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
//      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
//  
//  int err,status,ncid,groupid,cmode,val_id,r_id,retval;
//  cmode = NC_NOWRITE;
//#if USE_MPI==1
//  err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
//#else
//  err = nc_open(ncfilename.c_str(), cmode, &ncid);
//#endif
//  std::cout << "err " << err << std::endl;
//  std::cout << "Error: " <<  nc_strerror(err) << std::endl;
//  err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
//  std::cout << "Error: " <<  nc_strerror(err) << std::endl;
//  int ndimsp;
//  nc_inq_varndims (ncid, val_id, &ndimsp);
//  //Field::n_dimensions=ndimsp;
//  std::cout << "val number of dims " << ndimsp << std::endl;
//  int *dimidsp = new int[ndimsp];
//  nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
//  std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
//  std::vector<char> dimension_names(ndimsp); 
//  for(int i=0; i<ndimsp; i++)
//  {
//	  err = nc_inq_dimname(ncid,i, &dimension_names[i]);
//    std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
//  }
//
//  //size_t val_size=0;
//  //std::vector<float> coord1;
//  //std::vector<float> coord2;
//  //std::vector<float> coord3;
//
//  //for(int i=0;i<ndimsp;i++)
//  //{
//  //  size_t lenp;
//  //  nc_inq_dimlen(ncid,dimidsp[i],&lenp);
//  //  if(i==0)
//  //  {
//  //    val_size = lenp;
//  //  }
//  //  else
//  //  {
//  //    val_size = val_size*lenp;
//  //  }
//  //  std::cout << " i and dim size " << i << " " << lenp << std::endl;
//  //  char dimname[128];
//  //  nc_inq_dimname (ncid, dimidsp[i], &dimname[0]);
//  //  std::cout << " name " << dimname << std::endl;
//  //  int dim_var_id;
//  //  nc_inq_varid (ncid, dimname, &dim_var_id);
//  //  size_t startp=0;
//  //  //for (auto it = coord.cbegin(); it != coord.cend(); it++)
//  //  //{
//  //  //  std::cout << *it << ' ';
//  //  //}
//  //  if(i==0)
//  //  {
//  //    coord1.resize(lenp);
//  //    std::cout << "coord1 size " << coord1.size() << std::endl;
//  //    nc_get_vara_double (ncid, dim_var_id, &startp, &lenp, &coord1[0]);
//  //  for (auto it = coord1.cbegin(); it != coord1.cend(); it++)
//  //  {
//  //    std::cout << *it << ' ';
//  //  }
//  //  }
//  //  else if(i==1)
//  //  {
//  //    coord2.resize(lenp);
//  //    nc_get_vara_double (ncid, dim_var_id, &startp, &lenp, &coord2[0]);
//  //  }
//
//  //}
//  //std::vector<P> val(val_size);
//  //size_t startp=0;
//  //nc_get_var_double (ncid, val_id, &val[0]);
//  //std::cout << "val_id " << val_id  << std::endl;
//  //std::cout << "val_size " << startp << " " << val_size << std::endl;
//  ////for(int i=0 ;i<val_size; i++)
//  ////{
//  ////  std::cout << val[i] << " " ;
//  ////}
//  //nc_close(ncid);
//
//
//  const char* char_r = "r";
//  const char* char_z = "z";
//  if(ndimsp == 2 & dimension_names[0] == *char_r 
//		  & dimension_names[1] == *char_z)
//  {
//   FieldType this_field_type = FieldType_2d_xz;
//   field_type = this_field_type;
//   nD = 2;
//
//  }
  //  field_type{get_field_type(cfg, field_name)},	  
  //  nD{0}, dimensions{3,0},
  //  x{2,0.0},y{3,0.0},z{3,0.0},values{9,0.0},
  // fooHandler{returnPointerTable(1)} 
//    };
  //virtual float interpolate()
  //{
  //  std::cout << "field interpolate" << std::endl;
  //  return 0.0;
  //}
  CUDA_CALLABLE_MEMBER
  float returnOne(float x,float y, float z)
  {
	  return 1;
  }
  CUDA_CALLABLE_MEMBER
  float returnTwo(float x,float y, float z)
  {
	  float r = sqrt(x*x + y*y);
	  float value = interp_2d(r,z);
	  return value;
  }
  
  //CUDA_CALLABLE_MEMBER
  float interp_2d(float x,float y)
  {
	  float dx = this->x[1] - this->x[0];
	  float dy = this->y[1] - this->y[0];
          std::cout << "x y " << x << " " << y << std::endl;
          std::cout << "dx dy " << dx << " " << dy << std::endl;
	  int nx = this->dimensions[0];
	  int ny = this->dimensions[1];

	  int x_index = std::floor((x-this->x[0])/dx);
	  int y_index = std::floor((y-this->y[0])/dy);
          std::cout << "value 1 2 " << values[(x_index)*ny + y_index] << " " << values[(x_index+1)*ny + y_index] << std::endl;
          std::cout << "value 3 4 " << values[(x_index)*ny + y_index+1] << " " << values[(x_index+1)*ny + y_index+1] << std::endl;
          std::cout << "indices " << x_index << " " << y_index << std::endl;
	  float value_y0 = ((this->x[x_index + 1] - x)*this->values[x_index*ny + y_index] 
			  + (x - this->x[x_index])* this->values[(x_index + 1)*ny + y_index])/dx;
	  
	  float value_y1 = ((this->x[x_index + 1] - x)*this->values[x_index*ny + (y_index+1)] 
			  + (x - this->x[x_index])* this->values[(x_index + 1)*ny + (y_index+1)])/dx;
          std::cout << "value y0 y1 " << value_y0 << " " << value_y1 << std::endl;
	  float value = ((this->y[y_index + 1] - y)*value_y0 + (y - this->y[y_index] )*value_y1)/dy;
	  //value = this->values[x_index*ny + y_index];
	  //for (int i=0; i< this->x.size(); i++){
	  //        for(int j=0; j< this->y.size(); j++){
	  //      	  std::cout << "x " << this->x[i]
	  //      	            << "y " << this->y[j]
	  //      	            << "val " << this->values[i*ny+j]
	  //      		    << std::endl;
	  //        }
	  //}
	  return value;
  }
  
 // CUDA_CALLABLE_MEMBER
 // void BuildPointerTable(int num)
 //       {
 //       	if(num==0) this->fooHandler = &Field::returnOne;
 //       	else this->fooHandler = &Field::returnTwo;
 //       }
  
  CUDA_CALLABLE_MEMBER
  FunctionHandler<Field> returnPointerTable(int num)
	{
		if(num==FieldType_2d_xz){
			//std::cout << "foo1 " << std::endl;
		       	return &Field::returnOne;
		}
		else{
		       	return &Field::returnTwo;
			//std::cout << "foo2 " << std::endl;
		}
	}
  float interpolate(float x, float y, float z)
  {
	float val =   (this->*(this->fooHandler))(x,y,z);
		return val;
  }
  
  FieldType get_field_type(libconfig::Config &cfg,std::string field_name) 
  {
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      std::cout << "err " << err << std::endl;
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      //Field::n_dimensions=ndimsp;
      std::cout << "val number of dims " << ndimsp << std::endl;
      int dimidsp[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
      std::vector<char> dimension_names(ndimsp); 
      for(int i=0; i<ndimsp; i++)
      {
              err = nc_inq_dimname(ncid,i, &dimension_names[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
      }

      const char* char_r = "r";
      const char* char_z = "z";
      if(ndimsp == 2 & dimension_names[0] == *char_r 
            	  & dimension_names[1] == *char_z)
      {
       FieldType this_field_type = FieldType_2d_rz;
       return this_field_type;
      }
  }
  
  int get_field_n_dimensions(libconfig::Config &cfg,std::string field_name) 
  {
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      std::cout << "err " << err << std::endl;
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      //Field::n_dimensions=ndimsp;
      std::cout << "val number of dims " << ndimsp << std::endl;
      int dimidsp[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
      std::vector<char> dimension_names(ndimsp); 
      for(int i=0; i<ndimsp; i++)
      {
              err = nc_inq_dimname(ncid,i, &dimension_names[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
      }

      const char* char_r = "r";
      const char* char_z = "z";
      if(ndimsp == 2 & dimension_names[0] == *char_r 
            	  & dimension_names[1] == *char_z)
      {
       return 2;
      }
  }
  
  std::vector<std::size_t> get_field_dimensions(libconfig::Config &cfg,std::string field_name) 
  {
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      std::cout << "err " << err << std::endl;
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      //Field::n_dimensions=ndimsp;
      std::cout << "val number of dims " << ndimsp << std::endl;
      int dimidsp[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
      std::vector<char> dimension_names(ndimsp); 
      std::vector<std::size_t> dimension_lengths(ndimsp); 
      for(int i=0; i<ndimsp; i++)
      {
              err = nc_inq_dimname(ncid,i, &dimension_names[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
              err = nc_inq_dimlen(ncid,i, &dimension_lengths[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
      }

      const char* char_r = "r";
      const char* char_z = "z";
      if(ndimsp == 2 & dimension_names[0] == *char_r 
            	  & dimension_names[1] == *char_z)
      {
       return dimension_lengths;
      }
  }
  
  std::vector<std::string> get_field_dimension_names(libconfig::Config &cfg,std::string field_name) 
  {
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      std::cout << "err " << err << std::endl;
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      //Field::n_dimensions=ndimsp;
      std::cout << "val number of dims " << ndimsp << std::endl;
      int dimidsp[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
      std::vector<std::string> dimension_names(ndimsp); 
      std::vector<std::size_t> dimension_lengths(ndimsp); 
      for(int i=0; i<ndimsp; i++)
      {
	      char* name = new char[NC_MAX_NAME+1];
              err = nc_inq_dimname(ncid,i, name);
	      std::string name_string(name);
	      dimension_names[i] = name_string;
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
              err = nc_inq_dimlen(ncid,i, &dimension_lengths[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
      }

      const std::string char_r = "r";
      const std::string char_z = "z";
      if(ndimsp == 2 & dimension_names[0] == char_r 
            	  & dimension_names[1] == char_z)
      {
       return dimension_names;
      }
  }

  std::vector<float> get_field_vector(libconfig::Config &cfg,std::string field_name, int index) 
  {
	  std::cout << "get field vector " << index << std::endl;
      std::string ncfilename = get_variable<const char*>(cfg, field_name+".filename");
      std::string ncvarname = get_variable<const char*>(cfg, field_name+".variable_name");
  
      int err,status,ncid,groupid,cmode,val_id,r_id,retval;
      cmode = NC_NOWRITE;
#if USE_MPI==1
      err = nc_open_par(ncfilename.c_str(), cmode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
#else
      err = nc_open(ncfilename.c_str(), cmode, &ncid);
#endif
      std::cout << "err " << err << std::endl;
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      err=nc_inq_varid (ncid, ncvarname.c_str(), &val_id);
      std::cout << "Error: " <<  nc_strerror(err) << std::endl;
      int ndimsp;
      nc_inq_varndims (ncid, val_id, &ndimsp);
      //Field::n_dimensions=ndimsp;
      std::cout << "val number of dims " << ndimsp << std::endl;
      int dimidsp[ndimsp];
      nc_inq_vardimid(ncid,val_id, &dimidsp[0]);
      std::cout << "val_dim_ids " << dimidsp[0] << " " << dimidsp[1]<< std::endl;	  
      std::vector<std::string> dimension_names(ndimsp); 
      std::vector<std::size_t> dimension_lengths(ndimsp); 
      int val_size = 1;
      for(int i=0; i<ndimsp; i++)
      {
	      char* name = new char[NC_MAX_NAME+1];
              err = nc_inq_dimname(ncid,i, name);
	      std::string name_string(name);
	      dimension_names[i] = name_string;
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
              err = nc_inq_dimlen(ncid,i, &dimension_lengths[i]);
        std::cout << "dimension " << i << " name " << dimension_names[i] << std::endl;
	val_size = val_size*dimension_lengths[i];
      }
     std::cout << "dim lengths size " << dimension_lengths.size() << std::endl;
     if (index <= dimension_lengths.size() - 1 || index == -1 )
     {
       int vec_size;     
       std::string var_name;
       if( index == -1){
	       std::cout << "inside -1 " << std::endl;
	       var_name = ncvarname;
	       vec_size = val_size;
       }
       else if ( index == 0) {
	       var_name = dimension_names[index];   
	       vec_size = dimension_lengths[index];
       }    
       else if ( index == 1) {
	       var_name = dimension_names[index];   
	       vec_size = dimension_lengths[index];
       }    
       else if ( index == 2) {
	       var_name = dimension_names[index];   
	       vec_size = dimension_lengths[index];
       }    
       std::cout << "var name " << var_name << " " <<vec_size << std::endl; 
       int val_id;
       err=nc_inq_varid (ncid, var_name.c_str(), &val_id);
       std::cout << "var id " << val_id << std::endl; 
       std::vector<float> vec(vec_size,0.0);
       err = nc_get_var_float( ncid, val_id, &vec[0]);
       std::cout << "vec0 " << vec[2] << std::endl;
      nc_close(ncid);
       return vec;
       std::cout << "after return vec" << std::endl;
     }
     else{
	     std::vector<float> vec0(1,0.0);
	     return vec0;
     } 
  }

};

//class Field_2d_xz : public Field {
//public:
//  float interpolate()
//  {
//    std::cout << "field2dxz interpolate" << std::endl;
//    return 0.0;
//  };
//};
//class Field_2d_rz : public Field {
//public:
//  float interpolate()
//  {
//    std::cout << "field2drz interpolate" << std::endl;
//    return 0.0;
//  };
//};
//  Field* Field::Create(FieldType type) {
//    if (type == FieldType_2d_xz)
//        return new Field_2d_xz();
//    else if (type == FieldType_2d_xz)
//        return new Field_2d_xz();
//    else return NULL;
//}
  // Client class
 
class Field_client : public ManagedAllocation {
public:

    // Client doesn't explicitly create objects
    // but passes type to factory method "Create()"
    Field_client()
    {
        //FieldType type = FieldType_2d_xz;
	std::cout << "pfield"<< std::endl;
        auto pField = new Field();
	//std::cout << "pfield->create"<< std::endl;
	//pField->Create(type);
    }
    //Field* getField()  {
    //    return pField;
    //}

//private:
//    Field *pField;
};


#endif
