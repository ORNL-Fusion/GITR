gitrInput.cfg contains a *flags* struct. The values of the variables in this data field must
match their counterparts declared in GITR/CMake/user_options.cmake and set on the command line
during configuration with **-D**-style options. Example **flags** configuration in gitrInput.cfg:
```
flags = 
{
        USE_CUDA=1;
        USE_MPI=0;
        USEIONIZATION=0;
        USERECOMBINATION=0;
        USEPERPDIFFUSION=0;
        USECOULOMBCOLLISIONS=0;
        USETHERMALFORCE=0;
        USESURFACEMODEL=0;
        USESHEATHEFIELD=0;
        BIASED_SURFACE=0;
        USEPRESHEATHEFIELD=0;
        BFIELD_INTERP=0;
        LC_INTERP=0;
        GENERATE_LC=0;
        EFIELD_INTERP=0;
        PRESHEATH_INTERP=0;
        DENSITY_INTERP=0;
        TEMP_INTERP=0;
        FLOWV_INTERP=0;
        GRADT_INTERP=0;
        ODEINT=0;
        FIXEDSEEDS=0;
        PARTICLESEEDS=1;
        GEOM_TRACE=0;
        GEOM_HASH=0;
        GEOM_HASH_SHEATH=0;
        PARTICLE_TRACKS=1;
        PARTICLE_SOURCE_SPACE=0;
        PARTICLE_SOURCE_ENERGY=0;
        PARTICLE_SOURCE_ANGLE=0;
        PARTICLE_SOURCE_FILE=0;
        SPECTROSCOPY=0;
        USE3DTETGEOM=0;
        USECYLSYMM=0;
        USEFIELDALIGNEDVALUES=0;
        FLUX_EA=0;
        FORCE_EVAL=0;
        CHECK_COMPATIBILITY=1;
}
```
