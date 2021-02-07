import gitr
import io
import libconf
import numpy as np
from numpy import random
import netCDF4
import matplotlib.pyplot as plt

def setup_case(Zi = 1, amu_i = 2, Ti = 6,
               ni = 5e18, Bmag = 0.1, phi = 135,
               nP = 1e5, nT = 1e4, dt = 1.0e-9, height = 5.0e-4,
               filename = "input/gitrInput.cfg"):

    with io.open(filename) as f:
        config = libconf.load(f)
        config['backgroundPlasmaProfiles']['Z'] = Zi
        config['backgroundPlasmaProfiles']['amu'] = amu_i
        config['backgroundPlasmaProfiles']['Bfield']['r'] = Bmag*np.sin(np.deg2rad(phi))
        config['backgroundPlasmaProfiles']['Bfield']['z'] = Bmag*np.cos(np.deg2rad(phi))
        config['backgroundPlasmaProfiles']['Bfield']['y'] = 0.0
        config['backgroundPlasmaProfiles']['Temperature']['ti'] = Ti
        config['backgroundPlasmaProfiles']['Temperature']['te'] = Ti
        config['backgroundPlasmaProfiles']['Density']['ni'] = ni
        config['backgroundPlasmaProfiles']['Density']['ne'] = ni
        config['impurityParticleSource']['nP'] = int(nP)
        config['impurityParticleSource']['initialConditions']['impurity_Z'] = Zi
        config['impurityParticleSource']['initialConditions']['charge'] = Zi
        config['impurityParticleSource']['initialConditions']['impurity_amu'] = amu_i
        config['timeStep']['dt'] = dt
        config['timeStep']['nT'] = int(nT)

    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

    vTh = np.sqrt(2*Ti*1.602e-19/amu_i/1.66e-27);
    k = 1.38e-23*11604; 
    Bb = amu_i*1.66e-27/(2*Ti*k);
    vgrid = np.linspace(-3*vTh,3*vTh,1000);
    fv1 = np.sqrt(Bb/np.pi)*np.exp(-Bb*np.multiply(vgrid,vgrid));
    fv1CDF = np.cumsum(fv1);
    fv1CDF = np.divide(fv1CDF,fv1CDF[-1]);
    print('rands',random.rand(int(nP)))
    vx0 = np.interp(random.rand(int(nP)),fv1CDF,vgrid);
    vy0 = np.interp(random.rand(int(nP)),fv1CDF,vgrid);
    vz0 = np.abs(np.interp(random.rand(int(nP)),fv1CDF,vgrid));
    print('vx0',vx0)
    plt.hist(vx0, bins=30)
    plt.show()

    xdir = [np.cos(np.deg2rad(phi)), 0 , -np.sin(np.deg2rad(phi))];
    ydir = [0,1,0];
    zdir = [np.sin(np.deg2rad(phi)), 0 , np.cos(np.deg2rad(phi))];
    print('xdir', xdir)
    print('ydir', ydir)
    print('zdir', zdir)

    vx = xdir[0]*vx0 + ydir[0]*vy0 + zdir[0]*vz0;
    vy = xdir[1]*vx0 + ydir[1]*vy0 + zdir[1]*vz0;
    vz = xdir[2]*vx0 + ydir[2]*vy0 + zdir[2]*vz0;
    plt.hist(vz)
    plt.show()

    x = np.zeros(int(nP));
    y = np.zeros(int(nP));
    z = np.zeros(int(nP))+ height;

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

if __name__ == "__main__":
    setup_case(Zi = 1, amu_i = 1, Ti = 30,
               ni = 1e18, Bmag = 5, phi = 135,
               nP = 1e5, nT = 5e3, dt = 1.0e-9, height = 5.0e-4,
               filename = "input/gitrInput.cfg")
