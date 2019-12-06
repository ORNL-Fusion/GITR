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

#ifdef __CUDACC__
#include <curand_kernel.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/random.h>
#endif

#include <random>
enum FieldType { 
    FieldType_2d_xz,    FieldType_2d_rz 
};
class Field : public ManagedAllocation {
public:
  std::size_t nD;
  sim::Array<int> dimensions;
  sim::Array<float> x;
  sim::Array<float> y;
  sim::Array<float> z;
  sim::Array<float> values;
 static Field* Create(FieldType type); 
  Field() :
    nD{0}, dimensions{0,0},
    x{0,0.0},y{0,0.0},z{0,0.0},values{0,0.0} {};
  virtual float interpolate()
  {
    std::cout << "field interpolate" << std::endl;
    return 0.0;
  }
};

class Field_2d_xz : public Field {
public:
  float interpolate()
  {
    std::cout << "field2dxz interpolate" << std::endl;
    return 0.0;
  };
};
class Field_2d_rz : public Field {
public:
  float interpolate()
  {
    std::cout << "field2drz interpolate" << std::endl;
    return 0.0;
  };
};
  Field* Field::Create(FieldType type) {
    if (type == FieldType_2d_xz)
        return new Field_2d_xz();
    else if (type == FieldType_2d_xz)
        return new Field_2d_xz();
    else return NULL;
}
  // Client class
class Field_client {
public:

    // Client doesn't explicitly create objects
    // but passes type to factory method "Create()"
    Field_client()
    {
        FieldType type = FieldType_2d_xz;
        pField = Field::Create(type);
    }
    Field* getField()  {
        return pField;
    }

private:
    Field *pField;
};


#endif
