#python library tools for gitr

from distutils.dir_util import copy_tree
import sys
sys.path.append('/home/tqd/code/netcdf4-python')
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
def plot2dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)
    
    x1=np.array(config.geom.x1)
    x2=np.array(config.geom.x2)
    z1=np.array(config.geom.z1)
    z2=np.array(config.geom.z2)
    y1 = config.geom.y1
    y2 = config.geom.y2
    ys1 = np.ones(x1.size)*y1
    ys2 = np.ones(x1.size)*y2

    fig = plt.figure()
    #fig.patch.set_facecolor('black')
    ax = fig.gca(projection='3d')
    ax.plot(np.append(x1,x1[0]), np.append(ys1,ys1[0]), np.append(z1,z1[0]))
    ax.plot(np.append(x2,x2[0]), np.append(ys2,ys2[0]), np.append(z2,z2[0]))
    for i in range(0,x1.size):
        ax.plot(np.append(x1[i],x1[i]),np.append(y1,y2),np.append(z1[i],z1[i]))
    plt.savefig('geomPlot.png') 
    print('config', config) 
    print('x1 ', x1)
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
if __name__ == "__main__":
    nc_show("surface.nc")
    depositedEdist()
