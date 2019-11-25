#ifndef _FLAGS_
#define _FLAGS_
#include "libconfig.h++"
#include "utils.h"
class Flags
{
  private:
  public:
   int USE_IONIZATION;
   Flags();
   void initialize_flags(libconfig::Config &cfg);
};
#endif
