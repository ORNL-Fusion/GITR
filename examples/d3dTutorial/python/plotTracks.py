import netCDF4
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io,libconf

def plotTracks(filename='output/history.nc',geomFile='input/gitrGeometry.cfg',plot=1):
    
    ncFile = netCDF4.Dataset(filename,"r")
    nT = ncFile.dimensions['nT'].size
    nP = ncFile.dimensions['nP'].size
    x = np.reshape(np.array(ncFile.variables['x']),(nP,nT))
    y = np.reshape(np.array(ncFile.variables['y']),(nP,nT))
    z = np.reshape(np.array(ncFile.variables['z']),(nP,nT))
    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))
    
    with io.open(geomFile) as f:
        config = libconf.load(f)
    
    r_geom = config['geom']['x1']
    z_geom = config['geom']['z1']
    
    if plot ==1:
        plt.close() 
        plt.figure(1,figsize=(6, 10), dpi=60)
        if(x.shape[0] ==1):
            plt.plot(r[0][:],z[0][:],linewidth=5,color='green')
        else:
            for i in range(nP):
              plt.plot(r[i,:],z[i,:],linewidth=0.5)
        
        plt.autoscale(enable=True, axis='x', tight=True)
        r_geom = np.array(r_geom)
        z_geom = np.array(z_geom)
        plt.plot(np.append(r_geom,r_geom[0]),np.append(z_geom,z_geom[0]))
        plt.title('DIII-D W Impurity Simulation',fontsize=20)
        plt.xlabel('r [m]',fontsize=16)
        plt.ylabel('z [m]',fontsize=16)
        print('saving tracksRZ')
        plt.savefig('tracksRZ.png')
        plt.show()
if __name__ == "__main__":
    plotTracks(filename='output/history.nc',geomFile='input/gitrD3DGeometry2DWringsPy.cfg',plot=1)
