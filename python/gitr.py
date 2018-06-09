#python library tools for gitr

from distutils.dir_util import copy_tree
import sys
#sys.path.append('/home/tqd/code/netcdf4-python')
import netCDF4
import numpy as np
#import Tkinter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.use('agg')
#import cv2
import io,libconf
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pylab as pl
import scipy as sp
import math
import os
import shutil

def copy_folder(from_folder, to_folder):
    copy_tree(from_folder, to_folder)

def nc_show(filename):
    ncFile = netCDF4.Dataset(filename,"r")
    print("File : ",filename, " format ",ncFile.file_format)
    #print("nLines ", surfFile.dimensions.keys())
    print("DIMENSIONS: ")
    for key in ncFile.dimensions:
        print(key,ncFile.dimensions[key].size)

    print("VARIABLES: ")
    for key in ncFile.variables:
        print(key,ncFile.variables[key].dimensions,ncFile.variables[key][:])
def depositedEdist():
    ncFile = netCDF4.Dataset("surface.nc","r")
    nLines = ncFile.dimensions['nLines'].size
    nE = ncFile.dimensions['nAngles'].size
    nA = ncFile.dimensions['nEnergies'].size

    Z = np.array(ncFile.variables['Z'][:])
    Edist = np.array(ncFile.variables['surfEDist'])
    gridE = np.array(ncFile.variables['gridE'])
    gridA = np.array(ncFile.variables['gridA'])

    print('Z ', Z.shape)
    print('Edist ', Edist.shape)
    surfaceIndices = np.nonzero(Z)
    print('surfaces indices ', surfaceIndices)
    totalEdist = np.sum(Edist[surfaceIndices][:][:],axis=0)
    eadistSurf1 = Edist[0][:][:]
    edistSurf1 = np.sum(eadistSurf1, axis=0)
    print('total Edist',totalEdist.shape)
    EdistOnly = np.sum(totalEdist,axis=0)
    AdistOnly = np.sum(totalEdist,axis=1)
    totalImpacts = np.sum(EdistOnly,axis=0)
    print('total Impacts ', totalImpacts)
    E = np.linspace(0, 1000, 200)
    plt.figure(1,figsize=(10, 6), dpi=80)
    #plt.plot(E,edistSurf1)
    plt.plot(gridE,EdistOnly)
    #plt.show()
    plt.savefig('image1.png')
    plt.figure(2,figsize=(10, 6), dpi=80)
    plt.plot(gridA,AdistOnly)
    plt.savefig('image2.png')
    #image = cv2.imread("image1.png")
    #cv2.imshow("Image", image)
    for i in range(1,len(gridA)):
        writeEDfile("angle"+str(i)+".ED1",gridE,EdistOnly)

    writeEDfile("angularDist.dat",gridA,AdistOnly)
def writeEDfile(name,gridE,Edist):
    if(gridE[0] == 0):
        dE = gridE[1]-gridE[0]
        gridE = gridE+ 0.5*dE
    datafile_id = open(name, 'w+')
    data = np.array([gridE, Edist])
    data = data.T
    np.savetxt(datafile_id, data, fmt=['%.4f','%.4f'])
    #here the ascii file is populated. 
    datafile_id.close()
def modifyInputParam(filename="input/gitrInput.cfg",nT=100):
    with io.open(filename) as f:
        config = libconf.load(f)
    config['timeStep']['nT'] = nT    
    with io.open(filename,'w') as f:
        libconf.dump(config,f)
    
def plot2dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)
    
    x1=np.array(config.geom.x1)
    x2=np.array(config.geom.x2)
    z1=np.array(config.geom.z1)
    z2=np.array(config.geom.z2)
    Z=np.array(config.geom.Z)
    length=np.array(config.geom.length)
    y1 = config.geom.y1
    y2 = config.geom.y2
    ys1 = np.ones(x1.size)*y1
    ys2 = np.ones(x1.size)*y2
    #if plt.fignum_exists(num=1):
    plt.plot(np.append(x1,x1[0]),np.append(z1,z1[0]),linewidth=2.0,color='k')
    #else:    
    #    fig = plt.figure()
    #    #fig.patch.set_facecolor('black')
    #    ax = fig.gca(projection='3d')
    #    ax.plot(np.append(x1,x1[0]), np.append(ys1,ys1[0]), np.append(z1,z1[0]))
    #    ax.plot(np.append(x2,x2[0]), np.append(ys2,ys2[0]), np.append(z2,z2[0]))
    #    for i in range(0,x1.size):
    #        ax.plot(np.append(x1[i],x1[i]),np.append(y1,y2),np.append(z1[i],z1[i]))
    #    plt.savefig('geomPlot.png') 
    #    print('config', config) 
    #    print('x1 ', x1)
    return x1,x2,z1,z2,length,Z
def read3dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)
    
    x1=np.array(config.geom.x1)
    x2=np.array(config.geom.x2)
    x3=np.array(config.geom.x3)
    y1=np.array(config.geom.y1)
    y2=np.array(config.geom.y2)
    y3=np.array(config.geom.y3)
    z1=np.array(config.geom.z1)
    z2=np.array(config.geom.z2)
    z3=np.array(config.geom.z3)
    area=np.array(config.geom.area)
    Z = np.array(config.geom.Z)
    materialSurfaceInidces = np.nonzero(Z)
    surfIndArray = np.asarray(materialSurfaceInidces)
    print('Number of W surfaces ', surfIndArray.size)
    return x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray
def iter2dProcessing():
    x1,x2,z1,z2,length,Z = plot2dGeom('input/iter2dRefinedOuterTarget.cfg')
    plt.close()
    grossDep,grossEro,sumWeightStrike,E,A,EAdist,surfaceNumbers = nc_readSurface()
    plt.close()
    netErosion = grossDep - grossEro
    colormap = plt.cm.bwr
    normalize = matplotlib.colors.Normalize(vmin=-1.,vmax=1.) 
    plt.subplot(2,1,1)
    for i in range(0,len(surfaceNumbers)):
        plt.plot([x1[surfaceNumbers[i]], x2[surfaceNumbers[i]]], [z1[surfaceNumbers[i]], z2[surfaceNumbers[i]]],color=colormap(0.))
        plt.hold(True)
    plt.subplot(2,1,2)
    plt.scatter(surfaceNumbers,netErosion)
    plt.savefig('surfaces.png') 
    plt.close()
    x,y,r,z,charge = nc_plotPositions('output/positions.nc')
    plt.hist(charge)
    plt.savefig('charge.png')
    plt.close
    rSep = 5.5543
    zSep = -4.3970
    gridRsep = np.zeros(160)
    gridZsep = np.zeros(160)
    rMinusRsep = np.zeros(160)
    for i in range(160):
        print(i)
        thisSurf = surfaceNumbers[i+4]
        lastSurf = surfaceNumbers[i+3]
        print(thisSurf)
        print(z2[thisSurf])
        if(i == 0):
            rMinusRsep[i] = z2[thisSurf] - zSep
        else:
            rMinusRsep[i] = rMinusRsep[i-1] + length[lastSurf]
    plt.close('all')
    s1 = plt.scatter(rMinusRsep,grossEro[range(4,164)],c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep,grossDep[range(4,164)],c='red')
    s3 = plt.scatter(rMinusRsep,netErosion[range(4,164)],c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosion.png') 
    plt.close('all')
    s1 = plt.scatter(rMinusRsep,np.log10(grossEro[range(4,164)]),c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep,np.log10(grossDep[range(4,164)]),c='red')
    s3 = plt.scatter(rMinusRsep,np.sign(netErosion[range(4,164)])*np.log10(np.absolute(netErosion[range(4,164)])),c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosionlog.png') 
    plt.close()
def piscesProcessing(r=0.01,path=''):
    x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray = read3dGeom('input/gitrGeometryPisces1inch.cfg')
    r1 = np.sqrt(np.multiply(x1[surfIndArray],x1[surfIndArray]) + np.multiply(y1[surfIndArray],y1[surfIndArray]))
    r2 = np.sqrt(np.multiply(x2[surfIndArray],x2[surfIndArray]) + np.multiply(y2[surfIndArray],y2[surfIndArray]))
    r3 = np.sqrt(np.multiply(x3[surfIndArray],x3[surfIndArray]) + np.multiply(y3[surfIndArray],y3[surfIndArray]))
    grossDep,grossEro,sumWeightStrike,E,A,EAdist = nc_readSurface()
    condition = r1 < r
    dep = np.extract(condition,grossDep)
    ero = np.extract(condition,grossEro)
    strike = np.extract(condition,sumWeightStrike)
    areas = np.extract(condition,area)
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)

    backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
    backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
    time = float(config.postProcessing.time);
    nParticles = float(config.impurityParticleSource.nP);
    backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    erodedMass = time*backgroundIonsPerSec*184*1.66e-27*backgroundSputtYield*1000;
    erodedMassPerParticle = erodedMass/nParticles;
    netErosion = np.sum(ero - dep);
    netStrike = np.sum(strike)
    totalArea = np.sum(areas)
    impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec*backgroundSputtYield/nParticles;
    impurityFlux = netStrike/totalArea*impurityParticlePerSecondPerComputationalPartice;
    Wfrac = impurityFlux/backgroundFlux;
    Aweight = np.sum(EAdist,axis=0)
    print('W impurity flux ', impurityFlux)
    print('W impurity fraction ', Wfrac)
    #for i in surfIndArray:
    np.savetxt('gitrFluxE.dat', E)
    np.savetxt('gitrFluxAweight.dat', Aweight)
    np.savetxt('gitrFluxA.dat', A[:-1])
    np.savetxt('gitrFluxEAdist.dat', EAdist)
    
    if(path != ''): 
        np.savetxt(path+'/'+'gitrFluxE.dat', E)
        np.savetxt(path+'/'+'gitrFluxAweight.dat', Aweight)
        np.savetxt(path+'/'+'gitrFluxA.dat', A)
        np.savetxt(path+'/'+'gitrFluxEAdist.dat', EAdist)
    
    Dfrac = float(config.postProcessing.Dfrac);
    Hefrac = float(config.postProcessing.Hefrac);
    Tfrac = float(config.postProcessing.Tfrac);
    file = open('gitrOut.txt','w') 
    file.write('plasmaSpecies=He W D T\n') 
    file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
    file.write('flux='+str(backgroundFlux/1e18)+'\n') 
    file.write('gitrOutputDir='+os.getcwd()+'\n') 
    file.close() 

    if(path != ''): 
        shutil.copyfile('gitrOut.txt',path+'/'+'gitrOut.txt')
        #file = open(path+'/'+'gitrOut.txt','w') 
        #
        #file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n') 
        #file.write('flux='+str(backgroundFlux/1e18)+'\n') 
        #file.write('gitrOutputDir='+os.getcwd()+'\n') 
        #file.close() 

    E0=0.0;
    E=1000.0;
    nE=100;
    dE = E/nE;
    Evec = np.linspace(0.5*dE,E-0.5*dE,nE);
    A0=0.0;
    A=90.0;
    nA=90;
    for i in range(0,nA):
        print('Evec size ', Evec.size)
        print('EAdist[:,i] size ', EAdist[:,i].size)
        edOut = np.column_stack((Evec,EAdist[:,i]))
        np.savetxt('dist'+str(i)+'.dat', edOut)

def d3dProcessing(r=0.01,path=''):
    #plt.figure(1,figsize=(6, 10), dpi=1000)
    #plot2dGeom(filename='../input/gitrD3DGeometry2DWrings.cfg')
    grossDep,grossEro,sumWeightStrike,E,A,EAdist = nc_readSurface()
    print('W gross deposition ', grossDep)
    print('W gross erosion ', grossEro)
    #condition = r1 < r
    #dep = np.extract(condition,grossDep)
    #ero = np.extract(condition,grossEro)
    #strike = np.extract(condition,sumWeightStrike)
    #areas = np.extract(condition,area)
    #with io.open('input/gitrInput.cfg') as f:
    #    config = libconf.load(f)

    #backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
    #backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
    #time = float(config.postProcessing.time);
    #nParticles = float(config.impurityParticleSource.nP);
    #backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    #erodedMass = time*backgroundIonsPerSec*184*1.66e-27*backgroundSputtYield*1000;
    #erodedMassPerParticle = erodedMass/nParticles;
    #netErosion = np.sum(ero - dep);
    #netStrike = np.sum(strike)
    #totalArea = np.sum(areas)
    #impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec*backgroundSputtYield/nParticles;
    #impurityFlux = netStrike/totalArea*impurityParticlePerSecondPerComputationalPartice;
    #Wfrac = impurityFlux/backgroundFlux;
    #Aweight = np.sum(EAdist,axis=0)
    #print('W impurity flux ', impurityFlux)
    #print('W impurity fraction ', Wfrac)
def plot3dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)
    
    x1=np.array(config.geom.x1)
    x2=np.array(config.geom.x2)
    x3=np.array(config.geom.x3)
    y1=np.array(config.geom.y1)
    y2=np.array(config.geom.y2)
    y3=np.array(config.geom.y3)
    z1=np.array(config.geom.z1)
    z2=np.array(config.geom.z2)
    z3=np.array(config.geom.z3)
    xs=[]
    ys=[]
    zs=[]
    for i in range(0,x1.size-1):
        xs.append(x1[i])
        xs.append(x2[i])
        xs.append(x3[i])
        ys.append(y1[i])
        ys.append(y2[i])
        ys.append(y3[i])
        zs.append(z1[i])
        zs.append(z2[i])
        zs.append(z3[i])
    fig = plt.figure()
    #fig.patch.set_facecolor('black')
    #ax = fig.gca(projection='3d')
    #ax.plot_trisurf(xs,ys)
    print('xs ys zs', xs, ys, zs)
    ax = Axes3D(fig)
    verts = [zip(xs,ys,zs)]
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim3d(0.5,2.5)
    ax.set_ylim3d(-0.2,0.2)
    ax.set_zlim3d(-0.5,1.5)
    ax.add_collection3d(Poly3DCollection(verts))
    plt.savefig('geomPlot.png') 
def nc_plotHist(filename='history.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nT = ncFile.dimensions['nT'].size
    nP = ncFile.dimensions['nP'].size
    x = np.reshape(np.array(ncFile.variables['x']),(nP,nT))
    y = np.reshape(np.array(ncFile.variables['y']),(nP,nT))
    z = np.reshape(np.array(ncFile.variables['z']),(nP,nT))
    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))
#    print('z ', z[:,0])
    print('x ', x.shape)
    print('x ', x[0,:])
    print('r ', r.shape)
#    print('z ', z[0][:].shape)
#    #for i in range(nT):
#    #    print('r,z ',i, x[i][0],z[i][0]) 
#    single = x[0][:];
    plt.figure(1,figsize=(6, 10), dpi=60)
    #plot2dGeom(filename='../input/gitrGeometry.cfg')
    if(x.shape[0] ==1):
        plt.plot(r[0][:],z[0][:],linewidth=5,color='green')
    else:
        for i in range(nP):
          #print('i', i)  
#          print('size', r[:,i].size)  
#          print('r ', r[:,i])  
          plt.plot(r[i,:],z[i,:],linewidth=0.5)
#          #plt.plot(r[i,:],z[i,:],linewidth=1.0)
#          #plt.setp(linewidth=0.2)
    plt.xlim((5.3,6.0))
    plt.ylim((-4.4,-3.0))
    #plt.autoscale(enable=True, axis='x', tight=True)
    plt.title('DIII-D W Impurity Simulation',fontsize=20)
    plt.xlabel('r [m]',fontsize=16)
    plt.ylabel('z [m]',fontsize=16)
    #plt.axis('equal')
    print('saving tracksRZ')
    plt.savefig('tracksRZ.png')
    plt.show()
    plt.close()
    plt.figure(1,figsize=(10, 6), dpi=100)
    if(x.shape[0] ==1):
        plt.plot(x,y,linewidth=0.5)
    else:
        for i in range(80):
          #print('i', i)  
          plt.plot(x[i,:],y[i,:],linewidth=0.5)
    #plt.ylim((-0.02,0.02))
    #plt.xlim((-0.02,0.02))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.axis('equal')
    print('saving tracksXY')
    plt.savefig('tracksXY.png')
    return x,y,z,r
def nc_plotSpec(filename='spec.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nBins = ncFile.dimensions['nBins'].size
    nR = ncFile.dimensions['nR'].size
    nZ = ncFile.dimensions['nZ'].size
    n = np.array(ncFile.variables['n'])
    print('n ', n.shape)
    dens = n[nBins-1,:,:]
    print('dens ', dens.shape)
    plt.close()
    plt.figure(1,figsize=(10, 6), dpi=2000)
    plotsize = math.ceil(nBins**(0.5))
    for i in range(nBins-1,nBins):
        dens = np.log10(n[i,:,:])
        #plt.subplot(plotsize,plotsize,i+1)
        plot2dGeom('input/iter2dRefinedOuterTarget.cfg')
        plt.title("ITER W Impurity Density")
        plt.xlabel("r [m]")
        plt.ylabel("z [m]")
        plt.imshow(dens,extent =[4.0, 8.4, -4.6, 4.7],origin='lower')
        plt.colorbar(orientation='vertical')
    plt.savefig('image1.png')
def nc_plotSpec3D(filename='spec.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nBins = ncFile.dimensions['nBins'].size
    nR = ncFile.dimensions['nR'].size
    nY = ncFile.dimensions['nY'].size
    nZ = ncFile.dimensions['nZ'].size
    n = np.array(ncFile.variables['n'])
    print('n ', n.shape)
    dens = n[nBins-1,:,:,:]
    print('dens ', dens.shape)
    plt.close()
    plt.figure(1,figsize=(10, 6), dpi=2000)
    plt.imshow(n[0,5,:,:])
    plt.colorbar(orientation='vertical')
    #plotsize = math.ceil(nBins**(0.5))
    #for i in range(nBins):
    #    dens = np.log10(n[i,:,:])
    #    plt.subplot(plotsize,plotsize,i+1)
    #    plt.imshow(dens,origin='lower')
    #    plt.colorbar(orientation='vertical')
    plt.savefig('slice1.png')
def nc_plotPositions(filename='positions.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nP = ncFile.dimensions['nP'].size
    x = np.array(ncFile.variables['x'])
    y = np.array(ncFile.variables['y'])
    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))
    z = np.array(ncFile.variables['z'])
    charge = np.array(ncFile.variables['charge'])
    print('x ',x)
    print('y ',y)
    print('r ',r)
    print('z ',z)
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.scatter(x, y)
    ax = fig.add_subplot(122)
    ax.scatter(r,z)
    fig.savefig('positions.png')
    return x,y,r,z,charge
def nc_plotVz(filename='history.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nT = ncFile.dimensions['nT'].size
    nP = ncFile.dimensions['nP'].size
    x = np.reshape(np.array(ncFile.variables['x']),(nT,nP))
    y = np.reshape(np.array(ncFile.variables['y']),(nT,nP))
    z = np.reshape(np.array(ncFile.variables['z']),(nT,nP))
    vx = np.reshape(np.array(ncFile.variables['vx']),(nT,nP))
    vy = np.reshape(np.array(ncFile.variables['vy']),(nT,nP))
    vz = np.reshape(np.array(ncFile.variables['vz']),(nT,nP))
    vperp = np.sqrt(np.multiply(vx[nT-1,:],vx[nT-1,:]) + np.multiply(vy[nT-1,:],vy[nT-1,:]))
    pitchAngle = np.arctan(np.divide(vperp,vz[nT-1,:]))
    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))
    plt.figure(1,figsize=(10, 6), dpi=250)
    #plt.plot(z[:,1],vz[:,1],linewidth=0.5)
    plt.plot(z,vz,linewidth=0.5,label=str(np.linspace(0,nP-1,nP)))
    #plt.legend()
    #plt.savefig('vz.png')
    #for i in range(nP):
      #print('i', i) 
      #if((i > 20) and (i < 41)):
      #plt.plot(z[:,i],vz[:,i],linewidth=0.5,label=str(i*90.0/180))
      #print('size', r[:,i].size)  
##      print('r ', r[:,i])  
#      plt.plot(r[:,i],z[:,i],linewidth=0.5)
##      #plt.plot(r[i,:],z[i,:],linewidth=1.0)
##      #plt.setp(linewidth=0.2)
#    plt.savefig('tracksRZ.png')
#    plt.close()
#    plt.figure(1,figsize=(10, 6), dpi=1000)
#    for i in range(nP):
#      print('i', i)  
#      plt.plot(x[:,i],y[:,i],linewidth=0.5)
#    #plt.ylim((-5.0,5.0))
#    #plt.xlim((3.8,8.5))
#    plt.savefig('tracksXY.png')
    plt.title("10 eV Electron Phase Space Plot for Varying Pitch Angles")
    plt.xlabel("z [m]")
    plt.ylabel("v_par [m/s]")
    #plt.legend()
    plt.savefig('vz.png')
    plt.close()
    plt.figure(1,figsize=(10, 6), dpi=250)
    plt.hist(pitchAngle,bins=30)
    plt.savefig('pa.png')
def nc_readSurface(filename='output/surface.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nS = len(ncFile.dimensions['nSurfaces'])
    nE = len(ncFile.dimensions['nEnergies'])
    nA = len(ncFile.dimensions['nAngles'])
    E = np.linspace(0,1000.0,nE+1)
    A = np.linspace(0,90.0,nA+1)
    surfaceNumbers = np.array(ncFile.variables['surfaceNumber'])
    grossDep = np.array(ncFile.variables['grossDeposition'])
    grossEro = np.array(ncFile.variables['grossErosion'])
    sumWeightStrike = np.array(ncFile.variables['sumWeightStrike'])
    surfEdist = np.reshape(np.array(ncFile.variables['surfEDist']),(nS,nA,nE))
    surfEdist.reshape((nS,nA,nE))
    print('size surfEdist ', surfEdist.size)
    print('shape surfEdist ', surfEdist.shape)
    EAdist = np.sum(surfEdist,axis=0)
    print('shape EAdist ', EAdist.shape)
    EAdist = EAdist.reshape(nE,nA)
    plt.pcolor(EAdist)
    plt.colorbar()
    plt.savefig('EAdist.png')
    return grossDep,grossEro,sumWeightStrike,E,A,EAdist,surfaceNumbers
def plotPitch(filename='positions.nc'):
    ncFile = netCDF4.Dataset(filename,"r")
    nP = ncFile.dimensions['nP'].size
    x = np.array(ncFile.variables['x'])
    y = np.array(ncFile.variables['y'])
    z = np.array(ncFile.variables['z'])
    vx = np.array(ncFile.variables['vx'])
    vy = np.array(ncFile.variables['vy'])
    vz = np.array(ncFile.variables['vz'])
    vperp = np.sqrt(np.multiply(vx,vx) + np.multiply(vy,vy))
    pitchAngle = np.arctan(np.divide(vperp,vz))
    r = np.sqrt(np.multiply(x,x) + np.multiply(y,y))
    plt.figure(1,figsize=(10, 6), dpi=250)
    plt.hist(pitchAngle,bins=30)
    plt.savefig('pa.png')
    plt.close()
    plt.subplot(3,1,1)
    plt.hist(vx,bins=30)
    plt.subplot(3,1,2)
    plt.hist(vy,bins=30)
    plt.subplot(3,1,3)
    plt.hidast(vz,bins=30)
    plt.savefig('vs.png')
if __name__ == "__main__":
    #asdfanc_show("surface.nc")
    #depositedEdist()
    if(os.path.exists('output/history.nc')):
    	nc_plotHist('output/history.nc')
    if(os.path.exists('output/spec.nc')):
    	nc_plotSpec('output/spec.nc')
    iter2dProcessing()
    #nc_plotSpec3D()
    #nc_plotPositions()
    #nc_plotVz()
    #plotPitch()
    #piscesProcessing()
    #modifyInputParam()
    #nc_readSurface()
