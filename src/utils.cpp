#include "utils.h"
#include "libconfig.h++"
//#include "interp2d.hpp"
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
int make2dCDF(int nX, int nY, int nZ, float* distribution, float* cdf)
{
    int index=0;
    for(int i=0;i<nX;i++)
    {
      for(int j=0;j<nY;j++)
      {
        for(int k=0;k<nZ;k++)
        {
          index = i*nY*nZ + j*nZ + k;
          if(k==0)
          {
            cdf[index] = distribution[index];
          }
          else
          {
            cdf[index] = cdf[index-1] + distribution[index];
          }
        }  
      }  
    }
    for(int i=0;i<nX;i++)
    {
      for(int j=0;j<nY;j++)
      {
        if(cdf[i*nY*nZ + (j+1)*nZ - 1]>0.0)
        {
          for(int k=0;k<nZ;k++)
          {  
            index = i*nY*nZ + j*nZ + k;
            cdf[index] = cdf[index]/
                       cdf[index-k+nZ-1];
          }
        }
      }
    }
  return 0;
}
int regrid2dCDF(int nX, int nY, int nZ,float* xGrid,int nNew,float maxNew, float*cdf, float* cdf_regrid)
{
  //std::cout << " inside regrid function "<<nX << " " << nY << " " << nZ << std::endl;
  int lowInd=0;
  int index=0;
  float spline = 0.0;
  for(int i=0;i<nX;i++)
  {
    for(int j=0;j<nY;j++)
    {
      for(int k=0;k<nZ;k++)
      {
        index = i*nY*nZ + j*nZ + k;
        spline = interp1dUnstructured(xGrid[k],nNew,maxNew,&cdf[index-k],lowInd);
        if(std::isnan(spline) || std::isinf(spline)) spline = 0.0;
        cdf_regrid[index] = spline;  
        if(i==0 && j==0)
        {
          std::cout << "index xGrid[k] " << index << " " << xGrid[k] << " " << nNew << " " <<
              maxNew << " " << spline << std::endl;
        }
      }  
    }
  }
  return 0;  
}
