import io,libconf
import os
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import scipy.interpolate as scii

def create_particle_source(nP=1e3,filename='../input/profiles.nc',geom_file='gitrD3DGeometry2DWringsPy.cfg'):
    with io.open(geom_file) as f:
        config = libconf.load(f)
    
    r1_geom = np.array(config['geom']['x1'])
    z1_geom = np.array(config['geom']['z1'])
    r2_geom = np.array(config['geom']['x2'])
    z2_geom = np.array(config['geom']['z2'])
    surf_geom = np.array(config['geom']['surface'])

    surface_indices = np.where(surf_geom > 0) 
    surface_indices = surface_indices[0]
    
    print('r surface', r1_geom[surface_indices])
    print('z surface', z1_geom[surface_indices])
    nP = int(nP)
    x = 1.48406 + np.random.rand(1,int(nP))*(1.50302-1.48406)
    y = np.zeros(int(nP))
    z = 1.24722 - 0.5438*(x-1.48406)
    x = x-3e-3
    z = z-3e-3/0.5438
    vx = -3.0e3*np.ones(nP)*(0.5438/1.1383)
    vy = -3.0e3*np.ones(nP)
    vz = -3.0e3*np.ones(nP)*(1/1.1383)
    #x = []
    #z = []
    #for i in range(len(surf_ne)-1):
    #    nP_this_surf = int(nP*surf_ne[i]/ne_total)
    #    if(i == 0):
    #        nP_this_surf = nP_this_surf  + 12
    #    this_r = surf_midpoint_r[i] + (np.random.rand(1,nP_this_surf)-0.5)*(r1_geom[surface_indices[i+1]] - r1_geom[surface_indices[i]])
    #    print('nP and rs', nP_this_surf, this_r)
    #    x.extend(list(this_r[0]))
    #    z.extend(surf_midpoint_z[i]*np.ones(nP_this_surf))
    #print('x',len(x))    
    #print('z',len(z))    
    #x = np.array(x)    
    plt.scatter(x,z)
    plt.show()
    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    nr = rootgrp.createDimension("nP", nP)
    x_nc = rootgrp.createVariable("x","f8",("nP"))
    y_nc = rootgrp.createVariable("y","f8",("nP"))
    z_nc = rootgrp.createVariable("z","f8",("nP"))
    vx_nc = rootgrp.createVariable("vx","f8",("nP"))
    vy_nc = rootgrp.createVariable("vy","f8",("nP"))
    vz_nc = rootgrp.createVariable("vz","f8",("nP"))
    x_nc[:] = x
    y_nc[:] = 0.0*np.zeros(int(nP))
    z_nc[:] = z
    vx_nc[:] = vx
    vy_nc[:] = vy
    vz_nc[:] = 3.0e3*np.ones(int(nP))
    rootgrp.close()
if __name__ == "__main__":
    create_particle_source(nP=1e3,filename='../input/profiles.nc',geom_file='../input/gitrD3DGeometrySASPy.cfg')
