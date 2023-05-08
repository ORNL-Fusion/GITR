import numpy as np
import numpy.matlib
import netCDF4
import matplotlib.pyplot as plt
#generate particle source and profiles
def generate():
    ##################
    # INPUTS SECTION #
    ##################
    
    # Simulated particle mass (amu)
    m=184;
    
    # Domain half-length
    L=2.5; # meters
    
    # Flow speed toward the left boundary
    flowSpeed = -750;
    
    # Location at which plasma flow ends
    sv = 1.2; #meters
    
    # Ion temperature at the left boundary
    TD0 = 10;
    
    # Temperature gradient
    dT_ds = 1.38;
    
    # Number of simulated particles
    nP=1e6;
    
    # Plasma profile resolution
    nR = 10;
    nZ = 10001;
    
    #########################
    # END OF INPUTS SECTION #
    #########################
    
    ##################################
    # DEFINITIION OF PARTICLE SOURCE #
    ##################################
    
    # Initialize particle positions randomly between x0 to x1, z0 to z1, y=0
    x = 0.0*np.ones(int(nP));
    y = 0.0*np.ones(int(nP));
    z = 1.0e-6*np.ones(int(nP));
    
    # Initialize particle velocities to be thermal distribution at TD0
    nPoints = 200;
    maxE = 20;
    Eb = 8.79; # Surface binding energy
    a = 5; # Free parameter
    E = np.linspace(0,maxE,nPoints);
    dE = E[1];
    thompson2 = a*(a-1)*E*(Eb**(a-1))/(E+Eb)**(a+1);

    ecdf = np.cumsum(thompson2);
    ecdf = ecdf/ecdf[-1];

    
    # Sample 1D Maxwellian speed distribution in all velocity components
    rand1 = np.random.uniform(0,1,int(nP));
    
    randTheta = 2*np.pi*np.random.uniform(0,1,int(nP));
    randPhi = 0.5*np.pi*np.random.uniform(0,1,int(nP));
    Esamp = np.interp(rand1,ecdf,E);
    
    v = np.sqrt(2*Esamp*1.602e-19/m/1.66e-27);

    vx = v*np.cos(randTheta)*np.sin(randPhi);
    vy = v*np.sin(randTheta)*np.sin(randPhi);
    vz = v*np.cos(randPhi);
    
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
    
    ##################################
    # DEFINITIION OF PLASMA PROFILES #
    ##################################
    
    br = 2.25
    bz = 0.15
    angle = np.arctan(bz/br);
    print('Angle',angle)
    # Define 2D plasma profile grid
    r = np.linspace(-3000,3000,nR);
    z = np.linspace(0.0,2*L,nZ);
    
    # Flow velocity setup
    #vz = np.zeros(nZ);
    #vz[np.where(z<=sv)] = flowSpeed;
    #vz2D = np.matlib.repmat(vz,nR,1);
    
    # Ion temperature
    Ti =np.zeros(nZ);
    Ti = TD0+50*(-np.abs((z-L)**10) + L**10)/L**10
    Ti2D = np.matlib.repmat(Ti,nR,1);
    plt.close()
    plt.plot(z,Ti)
    
    plt.show()
    plt.savefig('temp.png')

    gradTi = bz/br*50/L**10*(-10*(z-L)**9)
    gradTi2D = np.matlib.repmat(gradTi,nR,1);
    plt.close()
    plt.plot(z,gradTi)
    
    plt.show()
    plt.savefig('gradTi.png')
   
    ni = 1e19*(6*(np.abs((z-L)**10))/L**10 + 1)
    ni2D = np.matlib.repmat(ni,nR,1);
    plt.close()
    plt.plot(z,ni)
    
    plt.show()
    plt.savefig('dens.png')

    v_para = 2e4*((z-L)**9)/2.5**9
    v2D = np.matlib.repmat(v_para,nR,1);
    vr2D = np.cos(angle)*v2D;
    vz2D = np.sin(angle)*v2D;
    plt.close()
    plt.plot(z,v_para)
    
    plt.show()
    plt.savefig('v_para.png')

    #vz2D = np.transpose(vz2D)
    Ti2D = np.transpose(Ti2D)
    ni2D = np.transpose(ni2D)
    vr2D = np.transpose(vr2D)
    vz2D = np.transpose(vz2D)
    gradTi2D = np.transpose(gradTi2D)
    
    # Write 2D plasma profiles to netcdf
    rootgrp = netCDF4.Dataset("input/profiles.nc", "w", format="NETCDF4")
    nrr = rootgrp.createDimension("nR", int(nR))
    nzz = rootgrp.createDimension("nZ", int(nZ))
    
    rrr = rootgrp.createVariable("gridR","f8",("nR"))
    zzz = rootgrp.createVariable("gridZ","f8",("nZ"))
    vzz = rootgrp.createVariable("vz","f8",("nZ","nR"))
    vyy = rootgrp.createVariable("vy","f8",("nZ","nR"))
    vxx = rootgrp.createVariable("vx","f8",("nZ","nR"))
    tii = rootgrp.createVariable("ti","f8",("nZ","nR"))
    gtix = rootgrp.createVariable("gradTiR","f8",("nZ","nR"))
    gtiy = rootgrp.createVariable("gradTiY","f8",("nZ","nR"))
    gtiz = rootgrp.createVariable("gradTiZ","f8",("nZ","nR"))
    gtex = rootgrp.createVariable("gradTeR","f8",("nZ","nR"))
    gtey = rootgrp.createVariable("gradTeY","f8",("nZ","nR"))
    gtez = rootgrp.createVariable("gradTeZ","f8",("nZ","nR"))
    tee = rootgrp.createVariable("te","f8",("nZ","nR"))
    nii = rootgrp.createVariable("ni","f8",("nZ","nR"))
    nee = rootgrp.createVariable("ne","f8",("nZ","nR"))
    rrr[:] = r
    zzz[:] = z
    vxx[:] = vr2D
    vyy[:] = 0*vr2D
    vzz[:] = vz2D
    tii[:] = Ti2D
    tee[:] = Ti2D
    nii[:] = ni2D
    nee[:] = ni2D
    gtix[:] = np.cos(angle)*(gradTi2D + 0)
    gtiy[:] = 0*gradTi2D
    gtiz[:] = np.sin(angle)*(gradTi2D + 0)
    rootgrp.close()
    
    ##########################
    # END OF PARTICLE SOURCE #
    ##########################
if __name__ == "__main__":
    generate()
