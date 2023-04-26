import netCDF4
import numpy as np
import matplotlib.pyplot as plt
#import pdb

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
flux_strength = 1;
dens_1d = flux_strength*dt/nP/dz*np.sum(n[nBins-1,:,:], axis=1)

plt.plot(z+0.5*dz,dens_1d)
plt.title("Case A C Impurity Density")
plt.xlabel("z [m]")
plt.ylabel("n [m-3]")
d0 = np.loadtxt("density_cuda.txt")
plt.plot(d0[:,1]+0.5*dz,d0[:,0])

plt.show()
plt.savefig('density.png')

np.savetxt("density.txt", np.c_[z+0.5*dz, dens_1d])
