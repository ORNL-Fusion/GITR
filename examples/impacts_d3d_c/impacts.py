import io
import libconf
import numpy as np
from numpy import random
import netCDF4
import matplotlib.pyplot as plt
import pandas as pd

def setup_case(Zi = 1, amu_i = 2, charge_i = 1,
               Zb = 1, amu_b = 2, charge_b = 1,
               Ti = 6,
               ni = 5e18, Te = 6, ne = 5e18, Bmag = 0.1, phi = 135,
               nP = 1e5, nT = 1e4, dt = 1.0e-9, height_factor = 4,
               filename = "input/gitrInput.cfg"):
    
    ME = 9.10938356e-31;
    MI = 1.6737236e-27;
    Q = 1.60217662e-19;
    EPS0 = 8.854187e-12;

    with io.open(filename) as f:
        config = libconf.load(f)
        config['backgroundPlasmaProfiles']['Z'] = charge_b
        config['backgroundPlasmaProfiles']['amu'] = amu_b
        config['backgroundPlasmaProfiles']['Bfield']['r'] = Bmag*np.sin(np.deg2rad(phi))
        config['backgroundPlasmaProfiles']['Bfield']['z'] = Bmag*np.cos(np.deg2rad(phi))
        config['backgroundPlasmaProfiles']['Bfield']['y'] = 0.0
        config['backgroundPlasmaProfiles']['Temperature']['ti'] = Ti
        config['backgroundPlasmaProfiles']['Temperature']['te'] = Te
        config['backgroundPlasmaProfiles']['Density']['ni'] = ni
        config['backgroundPlasmaProfiles']['Density']['ne'] = ne
        config['impurityParticleSource']['nP'] = int(nP)
        config['impurityParticleSource']['initialConditions']['impurity_Z'] = Zi
        config['impurityParticleSource']['initialConditions']['charge'] = charge_i
        config['impurityParticleSource']['initialConditions']['impurity_amu'] = amu_i
        config['timeStep']['dt'] = dt
        config['timeStep']['nT'] = int(nT)

    with io.open(filename, 'w') as f:
        libconf.dump(config, f)
    
    w = charge_i*Q*Bmag/(amu_i*MI);

    vTh = np.sqrt(2*Ti*1.602e-19/amu_i/1.66e-27);
    gyroradius = vTh/w
    print('gyroradius',gyroradius)

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
    #plt.hist(vx0, bins=30)
    #plt.show()

    xdir = [np.cos(np.deg2rad(phi)), 0 , -np.sin(np.deg2rad(phi))];
    ydir = [0,1,0];
    zdir = [np.sin(np.deg2rad(phi)), 0 , np.cos(np.deg2rad(phi))];
    print('xdir', xdir)
    print('ydir', ydir)
    print('zdir', zdir)

    vx = xdir[0]*vx0 + ydir[0]*vy0 + zdir[0]*vz0;
    vy = xdir[1]*vx0 + ydir[1]*vy0 + zdir[1]*vz0;
    vz = xdir[2]*vx0 + ydir[2]*vy0 + zdir[2]*vz0;
    #plt.hist(vz)
    #plt.show()

    x = np.zeros(int(nP));
    y = np.zeros(int(nP));
    z = np.zeros(int(nP))+ height_factor*gyroradius;

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
def WEST_case(location = 18, ion_charge = 1, filename = '/Users/tqd/Code/GITR/west/He_2d/WEST_He_OuterTarget.xlsx'):
    target = pd.read_excel(filename)
    print(target)
    print('ni',target['ni'][location])
    new_ni = 1.0e15*round(target['ni'][location]/1e15)
    new_ne = 1.0e15*round(target['ne'][location]/1e15)
    print('new_ni', new_ni)
    setup_case(Zi = 2, amu_i = 4, charge_i = ion_charge,
               Zb = 2, amu_b = 4, charge_b = target['ave_charge'][location],
               Ti = target['Ti'][location], ni = new_ni,
               Te = target['Te'][location], ne = new_ne,
               Bmag = target['b_mag'][location], phi = 180-target['b_angle'][location],
               nP = 1e5, nT = 5e3, dt = 1.0e-9,height_factor = 4,
               filename = "input/gitrInput.cfg")

def ProtoMPEx_case(surface_potential = 1000, T=4, n=5e18, sheath_factor = 1.0):
    setup_case(Zi = 1, amu_i = 2, charge_i = 1,
               Zb = 1, amu_b = 2, charge_b = 1,
               Ti = T, ni = n,
               Te = T, ne = n,
               Bmag = 0.1, phi = 180-85,
               nP = 1e4, nT = 5e5, dt = 5.0e-10,height_factor = 5*sheath_factor,
               filename = "input/gitrInput.cfg")
    
    with io.open('input/gitrGeometry.cfg') as f:
        config = libconf.load(f)
        config['geom']['potential'][0] = surface_potential

    with io.open('input/gitrGeometry.cfg', 'w') as f:
        libconf.dump(config, f)
    
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)
        config['backgroundPlasmaProfiles']['sheath_factor'] = sheath_factor
    
    with io.open('input/gitrInput.cfg', 'w') as f:
        libconf.dump(config, f)

def ProtoMPEx_case_O(surface_potential = 1000, T=4, n=5e18, sheath_factor = 1.0):
    setup_case(Zi = 8, amu_i = 16, charge_i = 1,
               Zb = 1, amu_b = 2, charge_b = 1,
               Ti = T, ni = n,
               Te = T, ne = n,
               Bmag = 0.1, phi = 180-85,
               nP = 1e4, nT = 2e6, dt = 1.0e-9,height_factor = 5*sheath_factor,
               filename = "input/gitrInput.cfg")
    
    with io.open('input/gitrGeometry.cfg') as f:
        config = libconf.load(f)
        config['geom']['potential'][0] = surface_potential

    with io.open('input/gitrGeometry.cfg', 'w') as f:
        libconf.dump(config, f)
    
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)
        config['backgroundPlasmaProfiles']['sheath_factor'] = sheath_factor
    
    with io.open('input/gitrInput.cfg', 'w') as f:
        libconf.dump(config, f)
def d3d_case_C(charge = 1, Tee=4, Tii = 3, n=5e18,btot = 1,br = 0.1, bt = 0.9, bz = -0.2, sheath_factor = 1.0):
    phi_angle = np.degrees(np.arccos(bz/btot));
    print('phi_angle',phi_angle)
    setup_case(Zi = 6, amu_i = 12, charge_i = charge,
               Zb = 1, amu_b = 2, charge_b = 1,
               Ti = Tii, ni = n,
               Te = Tee, ne = n,
               Bmag = 0.1, phi = 180-85,
               nP = 1e4, nT = 2e6, dt = 1.0e-9,height_factor = 5*sheath_factor,
               filename = "input/gitrInput.cfg")
    
if __name__ == "__main__":
    #WEST_case()
    #setup_case(Zi = 1, amu_i = 1, charge_i = 1,
    #           Zb = 1, amu_b = 1, charge_b = 1,
    #           Ti = 30, ni = 1e18, Te = 30, ne = 1e18,
    #           Bmag = 5, phi = 135,
    #           nP = 1e4, nT = 5e5, dt = 1.0e-11,
    #           height_factor = 5,
    #           filename = "input/gitrInput.cfg")
    #setup_case(Zi = 2, amu_i = 4, charge_i = 2,
    #           Zb = 2, amu_b = 4, charge_b = 1.5,
    #           Ti = 58.434091, ni = 2.6e18, Te = 23.889060, ne = 2.6e18,
    #           Bmag = 5.94, phi = 91.48,
    #           nP = 1e6, nT = 1e4, dt = 1.0e-9,
    #           height_factor = 10,
    #           filename = "input/gitrInput.cfg")
    ProtoMPEx_case()
