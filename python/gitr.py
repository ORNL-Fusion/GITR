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
    
if __name__ == "__main__":
    nc_show("surface.nc")
    depositedEdist()
