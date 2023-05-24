import netCDF4
import numpy as np
import matplotlib.pyplot as plt
#import pdb

def process_out():
    # Read spec file dimensions and variables
    filename='output/spec.nc'
    ncFile = netCDF4.Dataset(filename, "r")
    nBins = ncFile.dimensions['nBins'].size
    nR = ncFile.dimensions['nR'].size
    nZ = ncFile.dimensions['nZ'].size
    n = np.array(ncFile.variables['n'])
    nvz = np.array(ncFile.variables['nvz'])
    nE = np.array(ncFile.variables['nE'])
    z = np.array(ncFile.variables['gridZ'])
    nP = np.array(ncFile.variables['nP'])
    dt = np.array(ncFile.variables['dt'])
    
    # Calculate length
    dz = z[1] -z[0]
    
    dens = n[nBins - 1, :, :]
    plt.close()
    plt.figure(1, figsize=(10, 6), dpi=2000)
    flux_strength = 2*np.tanh(1);
    dens_1d = flux_strength*dt/nP/dz*np.sum(n[nBins-1,:,:], axis=1)
    flux_1d = flux_strength*dt/nP/dz*np.sum(nvz, axis=1)
    temp_1d = flux_strength*dt/nP/dz*np.sum(nE, axis=1)
    
    np.savetxt("dens_flux_temp.txt", np.c_[z+0.5*dz, dens_1d, flux_1d, temp_1d])

    return z, dens_1d, flux_1d, temp_1d
if __name__ == "__main__":
    process_out()
