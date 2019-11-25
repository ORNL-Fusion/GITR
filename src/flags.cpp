#include "flags.hpp"

Flags::Flags() : USE_IONIZATION(0) {}
void Flags::initialize_flags(libconfig::Config &cfg) 
{
  USE_IONIZATION = getVariable_cfg<int> (cfg,"flags.USE_IONIZATION");
}
