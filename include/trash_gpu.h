#include "array.h"

class operate_boris
{
  public:

  operate_boris( );

  CUDA_CALLABLE_MEMBER    
  void operator()(std::size_t indx); 

  double *data;
};

class run_boris
{
  public:

  run_boris();

  void run();

  operate_boris operate;

  sim::Array<double> array;
};
