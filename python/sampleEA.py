import netCDF4
import numpy as np
#import Tkinter
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.interpolate as scii
nParticles =1000
filename = 'angularDists.nc'
ncFile = netCDF4.Dataset(filename,"r")
nLoc = len(ncFile.dimensions['nLoc'])
nEgrid = len(ncFile.dimensions['nEgrid'])
nPhiGrid = len(ncFile.dimensions['nPhiGrid'])
nThetaGrid = len(ncFile.dimensions['nThetaGrid'])
Egrid = np.array(ncFile.variables['Egrid'])
phiGrid = np.array(ncFile.variables['phiGrid'])
thetaGrid = np.array(ncFile.variables['thetaGrid'])
Edist = np.reshape(np.array(ncFile.variables['Edist']),(nEgrid,nLoc))
phiDist = np.reshape(np.array(ncFile.variables['phiDist']),(nPhiGrid,nLoc))
thetaDist = np.reshape(np.array(ncFile.variables['thetaDist']),(nThetaGrid,nLoc))
plt.plot(Egrid,Edist)
plt.savefig('image1.png')
plt.plot(thetaGrid,thetaDist)
plt.savefig('image2.png')
plt.plot(phiGrid,phiDist)
plt.savefig('image3.png')
plt.close()
eCDF = np.cumsum(Edist,axis=0)
eCDF = eCDF/eCDF[-1]
phiCDF = np.cumsum(phiDist,axis=0)
phiCDF = phiCDF/phiCDF[-1]
thetaCDF = np.cumsum(thetaDist,axis=0)
thetaCDF = thetaCDF/thetaCDF[-1]
print('ecdf sshape',eCDF.shape)
plt.plot(Egrid,eCDF)
plt.savefig('image12.png')
plt.plot(thetaGrid,thetaCDF)
plt.savefig('image22.png')
plt.plot(phiGrid,phiCDF)
plt.savefig('image32.png')
plt.close()
rand2=np.random.rand(nParticles)
rand3=np.random.rand(nParticles)
rand4=np.random.rand(nParticles)
rmrsInd = 20
pEint = scii.interp1d(eCDF[:,rmrsInd],Egrid,bounds_error=False,fill_value=4.0)
pE = pEint(rand2)
pPhiint = scii.interp1d(phiCDF[:,rmrsInd],phiGrid,bounds_error=False,fill_value=0.0)
pPhi = pPhiint(rand3)
pThetaint = scii.interp1d(thetaCDF[:,rmrsInd],thetaGrid,bounds_error=False,fill_value=90.0)
pTheta = pThetaint(rand4)
plt.hist(pE)
plt.savefig('image4.png')
plt.close()
plt.hist(pPhi)
plt.savefig('image5.png')
plt.close()
plt.hist(pTheta)
plt.savefig('image6.png')
