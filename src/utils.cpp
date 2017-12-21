#include "utils.h"
#include "libconfig.h++"

using namespace std;
void checkFlags(libconfig::Config &cfg)
{
    std::cout << "Checking compatibility of compile flags with input file "
                      << std::endl;
    const char *flags0[] = {"flags.USE_CUDA","flags.USEMPI",
                            "flags.USE_BOOST","flags.USEIONIZATION",
                            "flags.USERECOMBINATION","flags.USEPERPDIFFUSION",
                            "flags.USECOULOMBCOLLISIONS",
                            "flags.USETHERMALFORCE","flags.USESURFACEMODEL",
                            "flags.USESHEATHEFIELD","flags.BIASED_SURFACE",
                            "flags.USEPRESHEATHEFIELD","flags.BFIELD_INTERP",
                            "flags.LC_INTERP","flags.GENERATE_LC", "flags.EFIELD_INTERP",
                            "flags.PRESHEATH_INTERP","flags.DENSITY_INTERP",
                            "flags.TEMP_INTERP",
                            "flags.FLOWV_INTERP","flags.GRADT_INTERP",
                            "flags.ODEINT","flags.FIXEDSEEDS",
                            "flags.PARTICLESEEDS","flags.GEOM_TRACE","flags.GEOM_HASH",
                            "flags.GEOM_HASH_SHEATH","flags.PARTICLE_TRACKS",
                            "flags.PARTICLE_SOURCE_SPACE",
                            "flags.PARTICLE_SOURCE_ENERGY",
                            "flags.PARTICLE_SOURCE_ANGLE",
                            "flags.SPECTROSCOPY","flags.USE3DTETGEOM","flags.USECYLSYMM",
                            "flags.FLUX_EA"};
        int flagValues[] =  {USE_CUDA, USEMPI, USE_BOOST,USEIONIZATION,
                             USERECOMBINATION,USEPERPDIFFUSION,USECOULOMBCOLLISIONS,
                             USETHERMALFORCE,USESURFACEMODEL,USESHEATHEFIELD,BIASED_SURFACE,
                             USEPRESHEATHEFIELD,BFIELD_INTERP,LC_INTERP, GENERATE_LC,
                             EFIELD_INTERP,
                             PRESHEATH_INTERP,DENSITY_INTERP,TEMP_INTERP,
                             FLOWV_INTERP,GRADT_INTERP,ODEINT,FIXEDSEEDS,
                             PARTICLESEEDS,GEOM_TRACE,GEOM_HASH,
                             GEOM_HASH_SHEATH,PARTICLE_TRACKS,PARTICLE_SOURCE_SPACE,
                             PARTICLE_SOURCE_ENERGY,PARTICLE_SOURCE_ANGLE,
                             SPECTROSCOPY,USE3DTETGEOM,USECYLSYMM,FLUX_EA};
            int check1;
            for (int i=0; i<sizeof(flagValues)/sizeof(int); i++)
               {
                  if(cfg.lookupValue(flags0[i], check1))
                  {
                     if (flagValues[i] != check1)
                     { std::cout << "incompatibility in " << flags0[i]
                                           << " between input file and binary" << std::endl;
                       exit(0);
                     }
                     else
                     {
                        std::cout << flags0[i] <<" = " << check1<< std::endl;
                     }
                  }
                  else
                  {
                     std::cout << flags0[i] <<" was not found" << std::endl;
                  }
              }
}

int getDimFromFile (libconfig::Config &cfg,const std::string& file,const std::string& section,
        const std::string& s)
{
  std::string str;
  getVariable(cfg,section+s,str);
  int dim = readFileDim(file,str);
  return dim;
}
