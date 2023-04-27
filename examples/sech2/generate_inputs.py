import numpy as np
import numpy.matlib
import netCDF4
#generate particle source and profiles
def generate():
    ##################
    # INPUTS SECTION #
    ##################
    
    # Simulated particle mass (amu)
    m=12;
    
    # Domain span
    z0=-1; # meters
    z1=1; # meters
    
    # Ion temperature at the left boundary
    TD = 10;
    
    # Number of simulated particles
    nP=1e6;
    
    #########################
    # END OF INPUTS SECTION #
    #########################
    
    ##################################
    # DEFINITIION OF PARTICLE SOURCE #
    ##################################
    
    # Initialize particle positions randomly between x0 to x1, z0 to z1, y=0
    z_grid = np.linspace(-1,1,10000);

    sech = np.divide(1,np.cosh(z_grid))
    pdf_sech2 = np.multiply(sech,sech);
    cdf = np.cumsum(pdf_sech2);
    cdf = cdf/cdf[-1];

    x = np.zeros(int(nP));
    y = np.zeros(int(nP));
    z = np.interp(np.random.uniform(0,1,int(nP)),cdf,z_grid);
    
    # Initialize particle velocities to be thermal distribution at TD
    # Create 1D Maxwellian spped distribution
    vTh = np.sqrt(2*TD*1.602e-19/m/1.66e-27);
    k = 1.38e-23*11604;
    B = m*1.66e-27/(2*TD*k);
    vgrid = np.linspace(-3*vTh,3*vTh,1000000);
    fv1 = np.sqrt(B/np.pi)*np.exp(-B*np.multiply(vgrid,vgrid));
    fv1CDF = np.cumsum(fv1);
    fv1CDF = fv1CDF/fv1CDF[-1];
    
    # Sample 1D Maxwellian speed distribution in all velocity components
    vx = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);
    vy = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);
    vz = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);
    
    # Write to netcdf file
    rootgrp = netCDF4.Dataset("input/particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", int(nP))
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = x
    yyy[:] = y
    zzz[:] = z
    vxx[:] = vx
    vyy[:] = vy
    vzz[:] = vz
    rootgrp.close()
    
    ##########################
    # END OF PARTICLE SOURCE #
    ##########################
    
if __name__ == "__main__":
    generate()
