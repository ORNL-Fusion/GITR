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
    
    plt.plot(z+0.5*dz,dens_1d)
    plt.title("Case A C Impurity Density")
    plt.xlabel("z [m]")
    plt.ylabel("n [m-3]")
    plt.show()
    plt.savefig('density.png')
    
    np.savetxt("density.txt", np.c_[z+0.5*dz, dens_1d])
if __name__ == "__main__":
    process_out()
