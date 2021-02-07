import tkinter
import matplotlib.pyplot as plt
import matplotlib
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    matplotlib.use('TkAgg')
elif platform == "darwin":
    # OS X
    matplotlib.use('TkAgg')
import numpy as np
import numpy.matlib as ml
import gitr
import scipy.interpolate as scii
import netCDF4
import os
from scipy.interpolate import griddata
import zipfile

def read_hpic_iead(num=1,file_path='/Users/tyounkin/Dissertation/ITER/mq3/final/hPIC_IEAD_solps_conditions/hPIC_IEAD_DATA'):
    iead = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_IEAD.dat', dtype='float',skiprows=0,delimiter=',')
    energy = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_EnergyCenters.dat', dtype='float',skiprows=0,delimiter=',')
    angle = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_AngleCenters.dat', dtype='float',skiprows=0,delimiter=',')
    print('Solps shape', iead.shape)
    print('energy', energy)
    print('angle', angle)

def plot_hpic_iead(solps_path='solpsTarg.txt',HpicDataFolder = '/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180730/hpicwork0004', \
         solpsIndex = [3]):
    os.mkdir('hpic_ieads')
    me   =  9.10938356e-31; # Electron mass
    mp   =  1.67262190e-27; # Proton mass
    
    # Specify the following three to generate the figures: 
    #   - File containing the SOLPS data (skip first row)
    #   - Folder containing hPIC data 
    #   - List of SOLPS points resolved with hPIC
    #   - Ai, Ion Atomic Mass
    #   - Zi, Charge State
    nSpecOut = len(solpsIndex)
    SOLPS = np.loadtxt(solps_path, dtype='float',skiprows=1,delimiter=',')
    solps_shape = SOLPS.shape
    print('solps shape',solps_shape[0],solps_shape[1])
    nL = solps_shape[0];
    nSolpsSpecies = int((solps_shape[1]-7)/2)
    solps_species = np.loadtxt('speciesList.txt', dtype='float',skiprows=1,delimiter=' ')
    Ai = solps_species[:,2]
    Zi = solps_species[:,3]
    SolpsLocationList = list(range(1,nL+1));
    #Ai  = [ 2, 3, 4, 4, 9, 9, 9, 9, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 ];
    #Zi  = [ 1, 1, 1, 2, 1, 2, 3, 4, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10 ];
    nS = len(Ai);
    #nA = 91;
    #nE = 91;
    
    for k in range(nSpecOut):
        energy = np.loadtxt(HpicDataFolder[k] +'/' +'SolpsPoint1' +'_EnergyCenters.dat', dtype='float',skiprows=0,delimiter=',')
        angle = np.loadtxt(HpicDataFolder[k] +'/' + 'SolpsPoint1' +'_AngleCenters.dat', dtype='float',skiprows=0,delimiter=',')
        nE = len(energy)
        nA = len(angle)
        IEADs = np.zeros([nL,nE,nA]);
        nn = np.zeros(nL);
        ff = np.zeros(nL);
        EEE = np.zeros([nL,nE]);
        AAA =  np.zeros([nL,nA]);
        for i in range(len(SolpsLocationList)):
            Ti  = SOLPS[i,3];
            Te  = SOLPS[i,4]; 
            fi  = SOLPS[i,5:(5+nSolpsSpecies)];
            ni  = SOLPS[i,5+nSolpsSpecies:5+2*nSolpsSpecies];
            print('ni', ni)
            nn[i] = ni[solpsIndex[k]];
            ff[i] = fi[solpsIndex[k]];
            ntot = sum(ni);
            ci = ni/ntot;
            print(Zi)
            print(ci)
            print(Ai)
            phi_floating = -Te*(np.log(sum(np.multiply(np.multiply(Zi,ci),np.sqrt(np.divide(2*np.pi*me/mp,Ai))))));
            phi_total = phi_floating + 0.5*Te;
            #EE   = np.linspace(0,24*Te,nE);
            #AA   = np.linspace(0,89,nA);
            #EEE[i,:] = EE;
            #AAA[i,:] = AA;
            ID   = 'SolpsPoint'+str(i+1);
            energy = np.loadtxt(HpicDataFolder[k] +'/' +ID +'_EnergyCenters.dat', dtype='float',skiprows=0,delimiter=',')
            angle = np.loadtxt(HpicDataFolder[k] +'/' +ID +'_AngleCenters.dat', dtype='float',skiprows=0,delimiter=',')
            EEE[i,:] = energy;
            AAA[i,:] = angle;
            IEAD = np.loadtxt( HpicDataFolder[k] +'/' +ID +'_IEAD.dat',dtype='float',delimiter = ',');
            IEADs[i,:,:] = IEAD
                #plt.close()
                #plt.pcolor(AA,EE,IEAD)
                #plt.colorbar(orientation='vertical')
                #plt.title('plot title')
                #plt.title('ITER Divertor Target Plasma Density')
                #plt.ylabel('n[m-3]')
                #plt.xlabel('R-Rsep[m]')
                #plt.savefig('iead00.png')
    
            rootgrp = netCDF4.Dataset("hpic_ieads/hpic"+str(k)+".nc", "w", format="NETCDF4")
            nc_const = rootgrp.createDimension("const", 1)
            nll = rootgrp.createDimension("nL", nL)
            #nss = rootgrp.createDimension("nS", nS)
            nee = rootgrp.createDimension("nE", nE)
            naa = rootgrp.createDimension("nA", nA)
            rr = rootgrp.createVariable("r","f8",("nL"))
            zzz = rootgrp.createVariable("z","f8",("nL"))
            aa = rootgrp.createVariable("A","f8",("const"))
            zz = rootgrp.createVariable("Z","f8",("const"))
            rmrs = rootgrp.createVariable("RmRs","f8",("nL"))
            gridee = rootgrp.createVariable("gridE","f8",("nL","nE"))
            gridaa = rootgrp.createVariable("gridA","f8",("nL","nA"))
            fff = rootgrp.createVariable("flux","f8",("nL"))
            nnn = rootgrp.createVariable("dens","f8",("nL"))
            ieadd = rootgrp.createVariable("iead","f8",("nL","nE","nA"))
            rr[:] = SOLPS[:,1]
            zzz[:] = SOLPS[:,2]
            aa[:] = Ai[solpsIndex[k]]
            zz[:] = Zi[solpsIndex[k]]
            rmrs[:] = SOLPS[:,0]
            gridaa[:] = AAA
            gridee[:] = EEE
            fff[:] = ff
            nnn[:] = nn
            ieadd[:] = IEADs
            rootgrp.close()
    
    zipf = zipfile.ZipFile('hpic_ieads.zip', 'w', zipfile.ZIP_DEFLATED)
    zipdir('hpic_ieads/', zipf)
    zipf.close()        
    #print('shapes', SOLPS[:,0].shape, ff.shape)
    #plt.close()
    #plt.plot(SOLPS[:,0],abs(ff.T))
    #plt.title('plot title')
    #plt.title('ITER Divertor Target Plasma Density')
    #plt.ylabel('n[m-3]')
    #plt.xlabel('R-Rsep[m]')
    #plt.yscale('log')
    #plt.savefig('flux.png')
    #plt.close()

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))

def plot_hpic_ieadDavide(solps_path='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180730/solpsTarg.txt',HpicDataFolder = '/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180730/hpicwork0004'):
    me   =  9.10938356e-31; # Electron mass
    mp   =  1.67262190e-27; # Proton mass
    
    # Specify the following three to generate the figures: 
    #   - File containing the SOLPS data (skip first row)
    #   - Folder containing hPIC data 
    #   - List of SOLPS points resolved with hPIC
    #   - Ai, Ion Atomic Mass
    #   - Zi, Charge State
    SOLPS = np.loadtxt(solps_path, dtype='float',skiprows=1,delimiter=',')
    nL = 36;
    SolpsLocationList = list(range(1,nL+1));
    Ai  = [ 2, 3, 4, 4, 9, 9, 9, 9, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 ];
    Zi  = [ 1, 1, 1, 2, 1, 2, 3, 4, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10 ];
    nS = len(Ai);
    nA = 90;
    nE = 240;
    nn = np.zeros([nS,nL]);
    ff = np.zeros([nS,nL]);
    EEE = np.zeros([nL,nE]);
    AAA =  np.zeros([nL,nA]);
    IEADs = np.zeros([nS,nL,nE,nA]);

    for i in range(len(SolpsLocationList)):
        Ti  = SOLPS[i,3];
        Te  = SOLPS[i,4]; 
        ni  = np.append(SOLPS[i,6], np.append(np.append(np.append(SOLPS[i,6], SOLPS[i,8:10]), SOLPS[i,11:15]), SOLPS[i,16:26]));
        print('ni', ni)
        fi  = SOLPS[i,[21+6,21+6,29,30,32,33,34,35,37,38,39,40,41,42,43,44,45,46]];
        nn[:,i] = ni;
        ff[:,i] = fi;
        ntot = sum(ni);
        ci = ni/ntot;
        phi_floating = -Te*(np.log(sum(np.multiply(np.multiply(Zi,ci),np.sqrt(np.divide(2*np.pi*me/mp,Ai))))));
        phi_total = phi_floating + 0.5*Te;
        EE   = np.linspace(0,24*Te,240);
        AA   = np.linspace(0,89.0,90);
        EEE[i,:] = EE;
        AAA[i,:] = AA;
        for j in range(len(Ai)):
            ID   = 'SolpsPoint'+str(i+1);
            print(('loop ID',ID))
            try:
                IEAD = np.loadtxt( HpicDataFolder +'/' +ID +'_IEAD_Sp'+ str(j) +'.dat',dtype='float',delimiter=',');
            except Exception as e:
                print(("Failed at preliminary EXCEPTION %s (check 'errors.txt')\n" % (str(e))))
                IEAD = np.loadtxt( HpicDataFolder +'/' +ID +'_IEAD_Sp'+ str(j) +'.dat',dtype='float');

            IEADs[j,i,:,:] = IEAD;
            #plt.close()
	    #plt.pcolor(AA,EE,IEAD)
            #plt.colorbar(orientation='vertical')
            #plt.title('plot title')
            #plt.title('ITER Divertor Target Plasma Density')
            #plt.ylabel('n[m-3]')
            #plt.xlabel('R-Rsep[m]')
            #plt.savefig('iead00.png')
    
    rootgrp = netCDF4.Dataset("hpic.nc", "w", format="NETCDF4")
    nll = rootgrp.createDimension("nL", nL)
    nss = rootgrp.createDimension("nS", nS)
    nee = rootgrp.createDimension("nE", nE)
    naa = rootgrp.createDimension("nA", nA)
    rr = rootgrp.createVariable("r","f8",("nL"))
    zzz = rootgrp.createVariable("z","f8",("nL"))
    aa = rootgrp.createVariable("A","f8",("nS"))
    zz = rootgrp.createVariable("Z","f8",("nS"))
    rmrs = rootgrp.createVariable("RmRs","f8",("nL"))
    gridee = rootgrp.createVariable("gridE","f8",("nL","nE"))
    gridaa = rootgrp.createVariable("gridA","f8",("nL","nA"))
    fff = rootgrp.createVariable("flux","f8",("nS","nL"))
    nnn = rootgrp.createVariable("dens","f8",("nS","nL"))
    ieadd = rootgrp.createVariable("iead","f8",("nS","nL","nE","nA"))
    rr[:] = SOLPS[:,1]
    zzz[:] = SOLPS[:,2]
    aa[:] = Ai
    zz[:] = Zi
    rmrs[:] = SOLPS[:,0]
    gridaa[:] = AAA
    gridee[:] = EEE
    fff[:] = ff
    nnn[:] = nn
    ieadd[:] = IEADs
    rootgrp.close()
    print('shapes', SOLPS[:,0].shape, ff.shape)
    plt.close()
    plt.plot(SOLPS[:,0],abs(ff.T))
    plt.title('plot title')
    plt.title('ITER Divertor Target Plasma Density')
    plt.ylabel('n[m-3]')
    plt.xlabel('R-Rsep[m]')
    plt.yscale('log')
    plt.savefig('flux.png')
    plt.close()

def readHpic(pathname = 'hpic_ieads',solps_inds=[3,4]):
    nSpec = len(solps_inds)
    nLocations = [];
    nE = []
    nA = []
    r = []
    z = []
    Z = []
    A = []
    RmRs = []
    gridE = []
    gridA = []
    flux = []
    dens = []
    IEAD = []

    for i in range(nSpec):
        fname = pathname + '/hpic'+str(i)+'.nc'

        ncFile = netCDF4.Dataset(fname,"r")
        nLocations.append(len(ncFile.dimensions['nL']))
        nE.append(len(ncFile.dimensions['nE']))
        nA.append(len(ncFile.dimensions['nA']))
        print("dimensions nL nE nA", nLocations, nE, nA)

        r.append(np.array(ncFile.variables['r']))
        z.append(np.array(ncFile.variables['z']))
        Z.append(np.array(ncFile.variables['Z']))
        A.append(np.array(ncFile.variables['A']))
        print("Z",Z)
        print("A",A)
        RmRs.append(np.array(ncFile.variables['RmRs']))
        gridE.append(np.array(ncFile.variables['gridE']))
        gridA.append(np.array(ncFile.variables['gridA']))
        print('gridE',gridE[i])
        print('gridA',gridA[i])
        #gridE = np.reshape(gridE,[nLocations,nE])
        #gridA = np.reshape(gridA,[nLocations,nA])
        flux.append(np.array(ncFile.variables['flux']))
        #flux =np.reshape(flux,[nSpecies,nLocations])
        dens.append(np.array(ncFile.variables['dens']))
        #dens =np.reshape(dens,[nSpecies,nLocations])
        
        IEAD.append(np.array(ncFile.variables['iead']))
        #IEAD = np.reshape(IEAD,[nSpecies,nLocations,nE,nA])
    print('iead shape',IEAD[0].shape)
    return RmRs,r,z,nE,nA,nLocations,nSpec,Z,A,dens,flux,gridE,gridA,IEAD
def readFtridynSelf(fname='ftridynSelf.nc'):
    ncFile = netCDF4.Dataset(fname,"r")
    nE = len(ncFile.dimensions['nE'])
    nA = len(ncFile.dimensions['nA'])
    nEdist = len(ncFile.dimensions['nEdistBins'])
    nAdist = len(ncFile.dimensions['nAdistBins'])
    print("dimensions nL nS nE nA", nE, nA, nEdist, nAdist)

    E = np.array(ncFile.variables['E'])
    A = np.array(ncFile.variables['A'])
    print("E",E)
    print("A",A)
    spyld = np.array(ncFile.variables['spyld'])
    spyld = np.reshape(spyld,[nE,nA])
    rfyld = np.array(ncFile.variables['rfyld'])
    rfyld = np.reshape(rfyld,[nE,nA])
    plt.pcolor(E,A,spyld.T)
    plt.colorbar(orientation='vertical')
    plt.title('plot title')
    plt.title('W on W Self Sputtering')
    plt.ylabel('A [degrees]')
    plt.xlabel('E[eV]')
    plt.savefig('WW.png')
    eDistEgrid = np.array(ncFile.variables['eDistEgrid'])
    cosDistAgrid = np.array(ncFile.variables['cosDistAgrid'])
    energyDist = np.array(ncFile.variables['energyDist'])
    energyDist = np.reshape(energyDist,[nEdist,nA,nE])
    plt.close()
    plt.plot(eDistEgrid,energyDist[-1,-3,:])
    plt.title('W on W Self Sputtering')
    plt.ylabel('A [degrees]')
    plt.xlabel('E[eV]')
    plt.savefig('Edist.png')
def readFtridynBackground(fname='ftridynBackground.nc'):
    ncFile = netCDF4.Dataset(fname,"r")
    nE = len(ncFile.dimensions['nE'])
    nA = len(ncFile.dimensions['nA'])
    nS = len(ncFile.dimensions['nS'])
    nEdist = len(ncFile.dimensions['nEdistBins'])
    nAdist = len(ncFile.dimensions['nAdistBins'])
    print("dimensions nL nS nE nA", nE, nA, nEdist, nAdist)

    E = np.array(ncFile.variables['E'])
    A = np.array(ncFile.variables['A'])
    print("E",E)
    print("A",A)
    spyld = np.array(ncFile.variables['spyld'])
    spyld = np.reshape(spyld,[nS,nE,nA])
    rfyld = np.array(ncFile.variables['rfyld'])
    rfyld = np.reshape(rfyld,[nS,nE,nA])
    for i in range(nS):
        plt.close()
        plt.pcolor(E,A,np.reshape(spyld[i,:,:],[nE,nA]).T)
        plt.colorbar(orientation='vertical')
        plt.title('plot title')
        plt.title('W on W Self Sputtering')
        plt.ylabel('A [degrees]')
        plt.xlabel('E[eV]')
        plt.savefig('Sputt'+str(i)+'.png')
    eDistEgrid = np.array(ncFile.variables['eDistEgrid'])
    phiGrid = np.array(ncFile.variables['phiGrid'])
    energyDist = np.array(ncFile.variables['energyDist'])
    energyDist = np.reshape(energyDist,[nS,nE,nA,nEdist])
    cosXdist = np.array(ncFile.variables['cosXDist'])
    cosXdist = np.reshape(cosXdist,[nS,nE,nA,nAdist])
    print(('Max edist and adist ', np.max(energyDist), np.max(cosXdist)))
    plt.close()
    plt.plot(eDistEgrid,energyDist[1,-1,0,:])
    plt.title('W on W Self Sputtering')
    plt.ylabel('A [degrees]')
    plt.xlabel('E[eV]')
    plt.savefig('Edist.png')
    return E,A,spyld,eDistEgrid,energyDist,phiGrid,cosXdist
def computeSputtYld(plots=False,hpic_zipfile='hpic_ieads',solps_inds = [3,4],nEdistPts=40,ftridyn_background_file = 'ftridynBackground.nc'):
    #plot_hpic_iead()
    nSpec = len(solps_inds)
    with zipfile.ZipFile(hpic_zipfile+'.zip', 'r') as zipObj:
        # Extract all the contents of zip file in current directory
        zipObj.extractall()
    Energy,Angle,spyld,eDistEgrid,energyDist,phiGrid,cosXdist=readFtridynBackground(ftridyn_background_file)
    RmRs,r,z,nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD=readHpic(hpic_zipfile,solps_inds)
        
    RmRs = RmRs[0]
    Y = np.zeros((2*nSpec,len(RmRs)))
    ee, aa = np.meshgrid(Energy, Angle)
    eDist_total = np.zeros((nSpec,len(RmRs),nEdistPts))
    aDist_total = np.zeros((nSpec,len(RmRs),len(phiGrid)))

    for i in range(nSpec):
        gridE_current = gridE[i]
        gridA_current = gridA[i]
        IEAD_current = IEAD[i]
        #print('size e a spyld', Energy.shape, Angle.shape, spyld.shape)
        #FIXME 1 index is for helium specifically
        eDist = energyDist[1,:,:,:]
        aDist = cosXdist[1,:,:,:]
        f = scii.interp2d(Energy, Angle, np.transpose(spyld[1,:,:]))
        if plots:
            plt.close()
            plt.pcolor(Energy, Angle, np.transpose(spyld[1,:,:]))
            plt.colorbar(orientation='vertical')
            plt.title('spyld')
            plt.ylabel('E[eV]')
            plt.xlabel('A[degrees]')
            plt.savefig('spyld.png')
            plt.close()
    

        flux_current =flux[i]
        Y[nSpec+i,:] = flux_current
    
        for j in range(len(RmRs)):
            print('For species', i, 'processing location',j)
            IEAD_current_location = IEAD_current[j,:,:]
            gridE_current_location = gridE_current[j,:]
            gridA_current_location = gridA_current[j,:]
            if plots:
                plt.close()
                plt.pcolor(gridE_current_location,gridA_current_location,np.transpose(IEAD_current_location))
                plt.colorbar(orientation='vertical')
                plt.title('IEAD')
                plt.ylabel('E[eV]')
                plt.xlabel('A[degrees]')
                plt.savefig('iead.png')
                plt.close()
            

            print('spyld shape',spyld.shape)
            spyld_current = f(gridE_current_location, gridA_current_location)
            spyld_current = gitr.interp_2d(Energy, Angle, spyld[1,:,:],gridE_current_location, gridA_current_location)
            #spyld_current = griddata((gridE_current_location.flatten(),gridA_current_location.flatten()),spyld.flatten(), (r_midpoint, z_midpoint), method='linear')
            
            if plots:
                plt.close()
                plt.pcolor(gridE_current_location,gridA_current_location, np.transpose(spyld_current))
                plt.colorbar(orientation='vertical')
                plt.title('Y')
                plt.ylabel('E[eV]')
                plt.xlabel('A[degrees]')
                plt.savefig('Ys.png')
                plt.close()

            Y_iead = np.multiply(IEAD_current_location,spyld_current)
            Y[i,j] = np.sum(Y_iead)/np.sum(IEAD_current_location)
     
            print('edist size',eDist.shape)
            
            if (Y[i,j] > 0.0):
                print('yield',Y[i,j])
                edist_current = gitr.interp_3d(Energy, Angle,eDistEgrid, eDist ,gridE_current_location, gridA_current_location,np.linspace(0.0,40.0,nEdistPts))
                phi_dist_current = gitr.interp_3d(Energy, Angle,phiGrid, aDist ,gridE_current_location, gridA_current_location,phiGrid)
                
                for ii in range(len(gridE_current_location)):
                    for jj in range(len(gridA_current_location)):
                        this_edist = edist_current[ii,jj,:]
                        edist_sum = np.sum(this_edist)
                        this_adist = phi_dist_current[ii,jj,:]
                        adist_sum = np.sum(this_adist)
                        if (edist_sum > 0.0):
                            this_edist = np.divide(this_edist,edist_sum)
                            eDist_total[i,j,:] = eDist_total[i,j,:] + Y_iead[ii,jj]*this_edist
                            this_adist = np.divide(this_adist,adist_sum)
                            aDist_total[i,j,:] = aDist_total[i,j,:] + Y_iead[ii,jj]*this_adist


            #print(gridE_current_location)
            #print(Energy[0])
            #erange = np.where(gridE_current_location > Energy[0])
            #erange_high = np.where(gridE_current_location > Energy[-1])
            #arange = np.where(gridA_current_location > Angle[-1])
            #gridA_current_location[arange] = Angle[-1]
            #gridE_current_location[erange_high] = Energy[-1]
            #print('erange',erange)
            #erange0 = erange[0]
            ##erange0 = erange0[0]
            #print('erange0',erange0)
            #if (len(erange0)):
            #    erange0 = erange0[0]
            #    print('erange0',erange0)
            #    print('Edist max', np.max(eDist))
            #    my_interpolating_function = scii.RegularGridInterpolator((np.log10(Energy),Angle,eDistEgrid),eDist,fill_value=0.0)
            #    pts = np.meshgrid(np.log10(gridE_current_location[erange0:]),gridA_current_location,eDistEgrid)
            #    pts_columns = np.zeros((len(gridE_current_location[erange0:])*len(gridA_current_location)*len(eDistEgrid),3))
            #    print('pts size',len(np.reshape(pts[0],-1)))
            #    print('pts_col size',pts_columns.shape)
            #    pts_columns[:,0] = np.reshape(pts[0],-1,order='F')
            #    pts_columns[:,1] = np.reshape(pts[1],-1,order='F')
            #    pts_columns[:,2] = np.reshape(pts[2],-1,order='F')
            #    #pts = np.array(pts)
            #    #pts = np.reshape(pts,(-1,3))
            #    print('pts size',pts_columns.shape)
            #    print(pts_columns)
            #    eDist_interp = my_interpolating_function(pts_columns)
            #    print('Edist interp', np.max(eDist_interp))
            #    eDist_interp = np.reshape(eDist_interp,(len(eDistEgrid),len(gridA_current_location),(len(gridE_current_location)-erange0)),order='F')
            #    eDist_full = np.zeros((len(eDistEgrid),len(gridA_current_location),len(gridE_current_location)))
            #    eDist_full[:,:,erange0:] = eDist_interp
            #    #eDist_full = np.zeros((len(gridE_current_location),len(gridA_current_location),len(eDistEgrid)))
            #    #eDist_full[erange0:,:,:] = eDist_interp
            #    plt.close()
            #    for ii in range(len(gridE_current_location)):
            #        for jj in range(len(gridA_current_location)):
            #            if ((gridE_current_location[ii] >= Energy[0]) and (gridE_current_location[ii] <= Energy[-1])):
            #                coord_array = np.zeros((len(eDistEgrid),3))
            #                coord_array[:,2] = eDistEgrid
            #                coord_array[:,1] = gridA_current_location[jj]*np.ones(len(eDistEgrid))
            #                coord_array[:,0] = np.log10(gridE_current_location[ii])*np.ones(len(eDistEgrid))
            #                plt.plot(my_interpolating_function(coord_array))
            #    
            #    plt.title('edists')
            #    plt.ylabel('E[eV]')
            #    plt.xlabel('A[degrees]')
            #    plt.savefig('edist.png')
            #    plt.close()

            #    eDist_weighted = np.zeros(len(eDistEgrid))
            #    print('Y_iead', Y_iead)
            #    print('iead sum',np.sum(IEAD_current_location))
            #    for ee in range(len(gridE_current_location)):
            #        for aa in range(len(gridA_current_location)):
            #            eDist_ee_aa_normalized = eDist_full[:,ee,aa]/np.sum(eDist_full[:,ee,aa])
            #            if (np.sum(eDist_full[:,ee,aa]) == 0):
            #                eDist_ee_aa_normalized = np.zeros(len(eDistEgrid))
            #            eDist_weighted = eDist_weighted + np.abs(flux_current[j])*Y_iead[ee,aa]*eDist_ee_aa_normalized/np.sum(IEAD_current_location)
            #    #eDist_full = np.zeros((len(gridE_current_location)*len(gridA_current_location)*len(eDistEgrid),3))
            #    #eDist_interp = scii.interpn([np.log10(Energy),Angle,eDistEgrid],eDist, [gridE_current_location,gridA_current_location,eDistEgrid])
            #else:
            #    eDist_weighted = np.zeros(len(eDistEgrid))

            #print('eDist weighted', eDist_weighted)
            #eDist_total[i,j,:] = eDist_weighted
    
    print('edist total',eDist_total[0,25,:])
    np.savetxt('Yields.txt', np.transpose(Y), delimiter=' ',header='Yields and fluxes')
    np.savetxt('Edists.txt', np.sum(eDist_total,axis=0), delimiter=' ',header='Edists')
    np.savetxt('Adists.txt', np.sum(aDist_total,axis=0), delimiter=' ',header='Adists')


def printBackgroundDist(path = '',rmrsPoints = [-0.1,0.02,0.09,0.2]):
    RmRs,r,zPoints,nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD=readHpic()
    flux = np.abs(flux)
    specNames = ['D','T','He','Ne']
    flux_fracs = np.zeros((len(specNames),len(rmrsPoints)))
    specInds = [0,1,[2,3],[8,9,10,11,12,13,14,15,16,17]]
    bgFlux = np.zeros((len(specNames),len(rmrsPoints)))
    for k in range(len(specNames)):
        if os.path.isdir("GITRoutput_"+specNames[k]):
            print("GITRoutput folder already exists")
        else:
            os.mkdir('GITRoutput_'+specNames[k])

    for i in range(len(rmrsPoints)):
        #idx = (np.abs(RmRs - rmrsPoints[i])).argmin()
        aa = RmRs - rmrsPoints[i]
        rsepPlus = min(ii for ii in aa if ii >= 0)
        print(("Position:", list(aa).index(rsepPlus)))
        print(("Value:", rsepPlus))
        idx = list(aa).index(rsepPlus)
        rsepPlus = RmRs[idx]
        rsepMinus = RmRs[idx-1]
        idxMinus = idx-1
        print(('index',idx))
        print(('RsepMinus and indxMInus',rsepMinus,idxMinus))
        thisDist = np.zeros((len(specNames),nE,nA))
        for j in range(len(specNames)):
            print(('flux shape',flux.shape))
            if (idx > 0) and (idx < len(RmRs)):
                bgFlux[j,i] = (rmrsPoints[i] - rsepMinus)/(rsepPlus-rsepMinus)*np.sum(flux[specInds[j],idx]) + (rsepPlus - rmrsPoints[i])/(rsepPlus-rsepMinus)*np.sum(flux[specInds[j],idxMinus])
            else:
                bgFlux[j,i] = 0.0
            if isinstance(specInds[j],int):
                if (idx > 0) and (idx < len(RmRs)):
                    thisDist[j,:,:] = (rmrsPoints[i] - rsepMinus)/(rsepPlus-rsepMinus)*IEAD[specInds[j],idx,:,:] + (rsepPlus - rmrsPoints[i])/(rsepPlus-rsepMinus)*IEAD[specInds[j],idxMinus,:,:]
            else:
                if (idx > 0) and (idx < len(RmRs)):

                    thisDist[j,:,:] = np.sum((rmrsPoints[i] - rsepMinus)/(rsepPlus-rsepMinus)*IEAD[specInds[j],idx,:,:] + (rsepPlus - rmrsPoints[i])/(rsepPlus-rsepMinus)*IEAD[specInds[j],idxMinus,:,:],axis=0)
                #thisDist[j,:,:] = np.sum(IEAD[specInds[j],idx,:,:],axis=0) #first dimension is angle, second is energy
        #thisDist = np.transpose(thisDist)
        print(('thisDistHe size', thisDist.shape))
        #print(thisDist)
        Aweight = np.sum(thisDist,axis=1)
        #gitrDir = 'GITRoutput/'+'gitr'+str(i)
        #if os.path.isdir(gitrDir):
        #    print("director already exists",gitrDir)
        #else:
        #    os.mkdir(gitrDir)
        
        for k in range(len(specNames)):
            gitrDirSpec = 'GITRoutput_'+specNames[k]+'/gitr'+str(i)
            if os.path.isdir(gitrDirSpec):
                print(("director already exists",gitrDirSpec))
            else:
                os.mkdir(gitrDirSpec)
    
            if isinstance(specInds[k],int):
                print('single charge species')
            else:		
                nChargeStates = len(specInds[k])
            for jj in range(nChargeStates):
                plt.close()
                plt.pcolor(gridA[idx,:],gridE[idx,:],IEAD[specInds[k][jj],idx,:,:])
                plt.colorbar(orientation='vertical')
                plt.title('IEAD'+str(jj+1))
                plt.ylabel('E[eV]')
                plt.xlabel('A[degrees]')
                plt.savefig(gitrDirSpec+'/iead'+str(jj+1)+'.png')
                plt.close()
                plt.pcolor(gridA[idx,:],gridE[idx,:],thisDist[k,:,:])
                plt.colorbar(orientation='vertical')
                plt.title('IEAD')
                plt.ylabel('E[eV]')
                plt.xlabel('A[degrees]')
                plt.savefig(gitrDirSpec+'/iead.png')
            thisGridE = gridE[idx,:]
            if thisGridE[0] == 0.0:
                thisGridE = thisGridE + 0.5*(thisGridE[1] - thisGridE[0])
            np.savetxt(gitrDirSpec+'/gitrFluxE.dat', thisGridE)
            np.savetxt(gitrDirSpec+'/gitrFluxAweight.dat', Aweight[k,:])
            np.savetxt(gitrDirSpec+'/gitrFluxA.dat',gridA[idx,:])
            np.savetxt(gitrDirSpec+'/gitrFluxEAdist.dat', thisDist[k,:,:])
            
            if(path != ''): 
                np.savetxt(gitrDirSpec+'/gitrFluxE.dat', thisGridE)
                np.savetxt(gitrDirSpec+'/gitrFluxAweight.dat', Aweight[k,:])
                np.savetxt(gitrDirSpec+'/gitrFluxA.dat',gridA[idx,:])
                np.savetxt(gitrDirSpec+'/gitrFluxEAdist.dat', thisDist[k,:,:])

            for j in range(0,nA):
                #print('Evec size ', gridE.size)
                #print('EAdist[:,i] size ', EAdist[:,i].size)
                edOut = np.column_stack((thisGridE,thisDist[k,:,j]))
                np.savetxt(gitrDirSpec+'/dist'+str(j)+'.dat', edOut)
    
        fluxTotal = np.sum(bgFlux[:,i])
        Dfrac =  bgFlux[0,i]/fluxTotal
        Hefrac = bgFlux[2,i]/fluxTotal
        Tfrac = bgFlux[1,i]/fluxTotal
        Nefrac = bgFlux[3,i]/fluxTotal
        if fluxTotal==0.0:
            Dfrac = 0.0
            Hefrac = 0.0
            Tfrac = 0.0
            Nefrac = 0.0
            Wfrac = 0.0
        flux_fracs[:,i] = [Hefrac,Dfrac,Tfrac,Nefrac]
        #file = open('gitrOut'+str(i)+'.txt','w') 
        #file.write('plasmaSpecies=He W D T Ne\n') 
        #file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+' '+str(Nefrac)+'\n') 
        #file.write('flux='+str(fluxTotal/1e18)+'\n') 
        #file.write('gitrOutputDir_W='+os.getcwd()+'gitrOutput_W'+'\n') 
        #file.write('gitrOutputDir_D='+os.getcwd()+'gitrOutput_D'+'\n') 
        #file.write('gitrOutputDir_T='+os.getcwd()+'gitrOutput_T'+'\n') 
        #file.write('gitrOutputDir_He='+os.getcwd()+'gitrOutput_He'+'\n') 
        #file.write('gitrOutputDir_Ne='+os.getcwd()+'gitrOutput_Ne'+'\n') 
        #file.close() 

        #if(path != ''): 
        #    shutil.copyfile('gitrOut.txt',path+'/'+'gitrOut.txt')
    return  bgFlux,flux_fracs
if __name__ == "__main__":   
    #plot_hpic_ieadDavide(solps_path='/project/projectdirs/atom/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180905/solpsTarg.txt',HpicDataFolder = '/project/projectdirs/atom/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180905/hpicdata0005')
    #plot_hpic_iead(solps_path='solpsTarg.txt', \
    #        HpicDataFolder = ['/Users/tyounkin/Dissertation/ITER/mq3/final/hPIC_IEAD_solps_conditions/hPIC_IEAD_DATA', \
    #        '/Users/tyounkin/Dissertation/ITER/mq3/final/hPIC_IEAD_He2_final/hPIC_IEAD_He2'],\
    #        solpsIndex = [3,4])
    #plot_hpic_ieadDavide()
    #nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD=readHpic()
    #readFtridynSelf()
    #readFtridynBackground()
    plot_hpic_iead()
    computeSputtYld()
    #printBackgroundDist()
