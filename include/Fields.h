#ifndef _FIELD_
#define _FIELD_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "array.h"
#include "managed_allocation.h"
#include "math.h"
#include <cstdlib>
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

class Field : public ManagedAllocation {
public:
  std::size_t nD;
  sim::Array<int> dimensions;
  sim::Array<float> x;
  sim::Array<float> y;
  sim::Array<float> z;
  sim::Array<float> values;
  
  Field() :
    nD{0}, dimensions{0,0},
    x{0,0.0},y{0,0.0},z{0,0.0},values{0,0.0} {};

  CUDA_CALLABLE_MEMBER
  Field(const std::size_t nD, std::string input_path, libconfig::Config &cfg,
        std::string cfg_string)
      :

        dimensions(getFieldNs(nD, input_path, cfg, cfg_string)),
        x{(nD == (std::size_t)0) ? 0 : dimensions[0], 0.0},
        y{(nD == (std::size_t)0) ? 0 : dimensions[1], 0.0},
        z{(nD == (std::size_t)0) ? 0 : dimensions[2], 0.0},
        values{(nD == (std::size_t)0)
                   ? getVariable(cfg, cfg_string + "x")
                   : getVarFromFile(cfg, "xCompString", cfg_string,
                                    dimensions[0] * dimensions[1] *
                                        dimensions[2])} {};

  CUDA_CALLABLE_MEMBER
  float interpolate();
  
  Field operator()(int a, int b, int c)
  {
        //Field F;
        this->nD=3;
	std::cout << "about to do something bad " << endl;
        this->dimensions.resize(nD),
	std::cout << "just did something bad " << endl;
        this->dimensions[0] = a;
	std::cout << "just did something else bad " << endl;
        this->dimensions[1] = b;
        this->dimensions[2] = c;
	std::cout << "done here " << endl;
	//return F;
  }
  CUDA_CALLABLE_MEMBER
  std::vector<int> getFieldNs(const std::size_t nD, std::string input_path,
                              libconfig::Config &cfg, std::string cfg_string) {

    int nx = 1, ny = 1, nz = 1;
    std::string file;
    importVectorFieldNs(cfg, input_path, nD, cfg_string, nx, ny, nz, file);
    std::vector<int> N(3);
    N[0] = nx;
    N[1] = ny;
    N[2] = nz;
    return N;
  }

  std::vector<float> getVariable(libconfig::Config &cfg, const std::string &s) {
    float val = 0.0;
    std::vector<float> value(1);
    if (cfg.lookupValue(s, val)) {
      // std::cout << s << " = " << tmp << std::endl;
    } else {
      std::cout << "ERROR: Failed importing " << s << std::endl;
      exit(0);
    }
    value[0] = val;
    return value;
  }

  std::string getVariableName(libconfig::Config &cfg, const std::string &s) {
    string val;
    if (cfg.lookupValue(s, val)) {
      std::cout << s << " = " << val << std::endl;
    } else {
      std::cout << "ERROR: Failed importing " << s << std::endl;
      exit(0);
    }
    return val;
  }

  std::vector<float> getVarFromFile(libconfig::Config &cfg,
                                    const std::string &variable,
                                    const std::string &section, int dim) {
    // Get NC file name from cfg
    std::string file = getVariableName(cfg, section + "fileString");
    // Get NC variable name from cfg
    std::string ncVarName = getVariableName(cfg, section + variable);
    // Get NC variable from NC file
    std::vector<float> values(dim);
    std::cout << "dim " << dim << std::endl;
    std::cout << "ncvarname " << ncVarName << std::endl;
    int dim2 = readFileVar(file, section, ncVarName, values[0]);
    return values;
  }

  //void broadcast(int world_rank)
  //{
  //  MPI_Bcast(&nR_Bfield,1,MPI_INT,0,MPI_COMM_WORLD);
  //  MPI_Bcast(&nY_Bfield,1,MPI_INT,0,MPI_COMM_WORLD);
  //  MPI_Bcast(&nZ_Bfield,1,MPI_INT,0,MPI_COMM_WORLD);
  //  MPI_Barrier(MPI_COMM_WORLD);
  //}
};

#endif
