import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as ml
import gitr
import scipy.interpolate as scii
import netCDF4
import os
from scipy.interpolate import griddata
def read_hpic_iead(num=1,file_path='/Users/tyounkin/Dissertation/ITER/mq3/final/hPIC_IEAD_solps_conditions/hPIC_IEAD_DATA'):
    iead = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_IEAD.dat', dtype='float',skiprows=0,delimiter=',')
    energy = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_EnergyCenters.dat', dtype='float',skiprows=0,delimiter=',')
    angle = np.loadtxt(file_path+'/SolpsPoint'+str(num)+'_AngleCenters.dat', dtype='float',skiprows=0,delimiter=',')
    print('Solps shape', iead.shape)
    print('energy', energy)
    print('angle', angle)

def plot_hpic_iead(solps_path='solpsTarg.txt',HpicDataFolder = '/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180730/hpicwork0004'):
    me   =  9.10938356e-31; # Electron mass
    mp   =  1.67262190e-27; # Proton mass
    
    # Specify the following three to generate the figures: 
    #   - File containing the SOLPS data (skip first row)
    #   - Folder containing hPIC data 
    #   - List of SOLPS points resolved with hPIC
    #   - Ai, Ion Atomic Mass
    #   - Zi, Charge State
    SOLPS = np.loadtxt(solps_path, dtype='float',skiprows=1,delimiter=' ')
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
        Te  = SOLPS[i,-4]; 
        ni  = np.append(SOLPS[i,5], np.append(np.append(np.append(SOLPS[i,5], SOLPS[i,7:9]), SOLPS[i,10:14]), SOLPS[i,15:25]));
        print('ni', ni)
        fluxInds = [28,30,32,33,35,36,37,38,40,41,42,43,44,45,46,47,48,49]
        fi  = SOLPS[i,fluxInds];
        nn[:,i] = ni;
        ff[:,i] = fi;
        ntot = sum(ni);
        ci = ni/ntot;
        phi_floating = -Te*(np.log(sum(np.multiply(np.multiply(Zi,ci),np.sqrt(np.divide(2*np.pi*me/mp,Ai))))));
        phi_total = phi_floating + 0.5*Te;
        EE   = np.linspace(0,24*Te,240);
        AA   = np.linspace(0,89,90);
        EEE[i,:] = EE;
        AAA[i,:] = AA;
        for j in range(len(Ai)):
            ID   = 'SolpsPoint'+str(i+1);
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

def readHpic(fname = 'hpic.nc'):
    ncFile = netCDF4.Dataset(fname,"r")
    nLocations = len(ncFile.dimensions['nL'])
    nSpecies = len(ncFile.dimensions['nS'])
    nE = len(ncFile.dimensions['nE'])
    nA = len(ncFile.dimensions['nA'])
    print("dimensions nL nS nE nA", nLocations, nSpecies, nE, nA)

    r = np.array(ncFile.variables['r'])
    z = np.array(ncFile.variables['z'])
    Z = np.array(ncFile.variables['Z'])
    A = np.array(ncFile.variables['A'])
    print("Z",Z)
    print("A",A)
    RmRs = np.array(ncFile.variables['RmRs'])
    gridE = np.array(ncFile.variables['gridE'])
    gridA = np.array(ncFile.variables['gridA'])
    gridE = np.reshape(gridE,[nLocations,nE])
    gridA = np.reshape(gridA,[nLocations,nA])
    flux = np.array(ncFile.variables['flux'])
    flux =np.reshape(flux,[nSpecies,nLocations])
    plt.close()
    plt.plot(RmRs,abs(flux.T))
    plt.title('plot title')
    plt.title('ITER Divertor Target Plasma Density')
    plt.ylabel('n[m-3]')
    plt.xlabel('R-Rsep[m]')
    plt.yscale('log')
    plt.savefig('flux.png')
    plt.close()
    dens = np.array(ncFile.variables['dens'])
    dens =np.reshape(dens,[nSpecies,nLocations])
    
    IEAD = np.array(ncFile.variables['iead'])
    IEAD = np.reshape(IEAD,[nSpecies,nLocations,nE,nA])
    print("gridA ", gridA[0,:]) 
    print("gridE ", gridE[0,:]) 
    print("iead shape ", IEAD.shape)
    print("iead00 ", np.reshape(IEAD[0,0,:,:],(nE,nA)))
    plt.pcolor(gridA[0,:],gridE[0,:],np.reshape(IEAD[0,0,:,:],(nE,nA)))
    plt.colorbar(orientation='vertical')
    plt.title('plot title')
    plt.title('ITER Divertor Target Plasma Density')
    plt.ylabel('n[m-3]')
    plt.xlabel('R-Rsep[m]')
    plt.savefig('iead00.png')
    #for i in range(len(z)):
    #    idx = (np.abs(zPoints - z[i])).argmin()
    #    print('index',idx)
    #    heFlux[i] = flux[idx] 
    #    thisDist = EAdist[:,:,idx] #first dimension is angle, second is energy
    #    thisDist = np.transpose(thisDist)
    #    print('thisDistHe size', thisDist.shape)
    #    #print(thisDist)
    #    Aweight = np.sum(thisDist,axis=0)
    #    gitrDir = 'GITRoutput/'+'gitr'+str(i)
    #    os.mkdir(gitrDir)
    #    np.savetxt(gitrDir+'/gitrFluxE.dat', gridE)
    #    np.savetxt(gitrDir+'/gitrFluxAweight.dat', Aweight)
    #    np.savetxt(gitrDir+'/gitrFluxA.dat',gridA[:-1])
    #    np.savetxt(gitrDir+'/gitrFluxEAdist.dat', thisDist)
    #    
    #    if(path != ''): 
    #        np.savetxt(path+'/'+gitrDir+'/gitrFluxE.dat', gridE)
    #        np.savetxt(path+'/'+gitrDir+'/gitrFluxAweight.dat', Aweight)
    #        np.savetxt(path+'/'+gitrDir+'/gitrFluxA.dat', gridA)
    #        np.savetxt(path+'/'+gitrDir+'/gitrFluxEAdist.dat', thisDist)

    #    for j in range(0,nA):
    #        #print('Evec size ', gridE.size)
    #        #print('EAdist[:,i] size ', EAdist[:,i].size)
    #        edOut = np.column_stack((gridE,thisDist[:,j]))
    #        np.savetxt(gitrDir+'/dist'+str(j)+'.dat', edOut)
    return RmRs,r,z,nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD
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
    energyDist = np.reshape(energyDist,[nE,nA,nEdist])
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
    #plt.close()
    #plt.plot(eDistEgrid,energyDist[-1,-3,:])
    #plt.title('W on W Self Sputtering')
    #plt.ylabel('A [degrees]')
    #plt.xlabel('E[eV]')
    #plt.savefig('Edist.png')
    return E,A,spyld,eDistEgrid,energyDist,phiGrid,cosXdist
def computeSputtYld(plots=0):
    #plot_hpic_iead()
    Energy,Angle,spyld,eDistEgrid,energyDist,phiGrid,cosXdist=readFtridynBackground()
    RmRs,r,z,nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD=readHpic()
    ftSpecies=0
    specInd = []
    yfElements = np.zeros((4,nLocations,nE,nA))
    yfSpecies = 0*IEAD
    sputtContribution = np.zeros([nLocations,len(A)])
    fContribution = np.zeros([nLocations,len(A)])
    for i in range(nLocations):
        if(gridE[i,-1] > 100.0):
            ftSpecies=0
            for j in range(len(A)):
                if plots:
                    plt.close()
                    plt.figure(1,figsize=(10, 6), dpi=200)
                    plt.subplot(2,2,1)
                    plt.pcolor(gridA[i,:],gridE[i,:],np.reshape(IEAD[j,i,:,:],(nE,nA)))
                    plt.colorbar(orientation='vertical')
                    plt.subplot(2,2,2)
                    plt.pcolor(Angle,Energy,np.reshape(spyld[ftSpecies,:,:],[len(Energy),len(Angle)]))
                    plt.colorbar(orientation='vertical')
                print('doing meshgrid')
                ieadGridX,ieadGridY = np.meshgrid(gridE[i,:],gridA[i,:])
                points = np.zeros([len(Energy)*len(Angle),2])
                val1 = np.reshape(spyld[ftSpecies,:,:],[len(Energy),len(Angle)])
                values = []
                print('doing points and values')
                for ii in range(len(Energy)):
                    for jj in range(len(Angle)):
                        points[ii*len(Angle)+jj,0] = Energy[ii]
                        points[ii*len(Angle)+jj,1] = Angle[jj]
                        values.append(val1[ii,jj])
                print('doing griddatta')
                grid_z1 = griddata(points, values, (ieadGridX, ieadGridY), method='linear',fill_value=0.0)
                print('grid_z',grid_z1.shape)
                if plots:
                    plt.subplot(2,2,3)
                    print('plotting gridz',gridE)
                    plt.pcolor(ieadGridY,ieadGridX,grid_z1)
                    plt.colorbar(orientation='vertical')
                print('saving')
                mult = np.multiply(grid_z1,np.reshape(IEAD[j,i,:,:],(nE,nA)).T)
                yfSpecies[j,i,:,:] = mult.T
                if plots:
                    plt.subplot(2,2,4)
                    plt.pcolor(ieadGridY,ieadGridX,mult)
                    plt.colorbar(orientation='vertical')
                    plt.savefig('hpicFT'+str(i)+str(j)+'.png')
                if(j < len(A)-1):
                    if(ftSpecies ==0):
                        specInd.append(j)
                        ftSpecies = ftSpecies+1
                    elif(A[j+1] > A[j]):
                        specInd.append(j)
                        print("aj+1 aj ftSpecies", A[j+1] ,A[j],ftSpecies)
                        ftSpecies = ftSpecies+1
                sputtContribution[i,j] = np.sum(mult)
                fContribution[i,j] = np.sum(np.reshape(IEAD[j,i,:,:],(nE,nA)))
    normSputtCont = np.divide(sputtContribution,np.sum(fContribution,axis=1)[:,None])
    normSputtCont[np.isnan(normSputtCont)] = 0.0
    yfElements[0,:,:,:] = yfSpecies[0,:,:,:]
    yfElements[1,:,:,:] = yfSpecies[1,:,:,:]
    yfElements[2,:,:,:] = np.sum(yfSpecies[2:4,:,:,:],axis=0)
    yfElements[3,:,:,:] = np.sum(yfSpecies[8:18,:,:,:],axis=0)
    print(('max yfSpc and yfEl', np.max(yfSpecies[8:18,:,:,:]), np.max(yfElements[3,:,:,:])))
    eDistLocations = np.zeros((nLocations,len(eDistEgrid)))
    aDistLocations = np.zeros((nLocations,len(phiGrid)))
    for i in range(nLocations):
        for j in range(nE):
            thisE = gridE[i,j]
            aa = Energy - thisE
            indPlusVal = min(ii for ii in aa if ii >= 0)
            idxPlus = list(aa).index(indPlusVal)
            eInd = idxPlus-1
            #print('Energy and eInd',thisE,eInd)
            if eInd >= 0:
                for k in range(nA):
                    thisA = gridA[i,k]
                    #print('Energy angle',thisE,thisA)
                    aInd = int((thisA - Angle[0])/(Angle[1]-Angle[0]))
                    #print('Energy angle ind',eInd,aInd)
                    #if yfElements[3,i,j,k] > 0.0:
                        #print('yfElements and energyDist',yfElements[3,i,j,k],energyDist[3,eInd,aInd,:])
                    eDistLocations[i,:] = yfElements[3,i,j,k]*energyDist[3,eInd,aInd,:] + eDistLocations[i,:]
                    aDistLocations[i,:] = yfElements[3,i,j,k]*cosXdist[3,eInd,aInd,:] + aDistLocations[i,:]
    rootgrp = netCDF4.Dataset("EAdistLocations.nc", "w", format="NETCDF4")
    nll = rootgrp.createDimension("nL", nLocations)
    nee = rootgrp.createDimension("nEdist", len(eDistEgrid))
    naa = rootgrp.createDimension("nAdist", len(phiGrid))
    eee = rootgrp.createVariable("E","f8",("nEdist"))
    aaa = rootgrp.createVariable("phi","f8",("nAdist"))
    edl = rootgrp.createVariable("eDistLocations","f8",("nL","nEdist"))
    adl = rootgrp.createVariable("aDistLocations","f8",("nL","nAdist"))
    eee[:] = eDistEgrid
    aaa[:] = phiGrid
    edl[:] = eDistLocations
    adl[:] = aDistLocations
    rootgrp.close()
                    
    plt.close()
    #plt.plot(RmRs,normSputtCont[:,0:8])
    plt.plot(RmRs,normSputtCont[:,0:4])
    plt.plot(RmRs,normSputtCont[:,8:16],linestyle='dashed')
    plt.plot(RmRs,normSputtCont[:,16:19],linestyle='dotted')
    plt.ticklabel_format(style='sci', axis='y',scilimits=(-4,-4))
    plt.title('Energy-Angle Averaged Sputtering Yield \n of Background Species\n')
    plt.ylabel('Yield')
    plt.xlabel('R-Rsep [m]')
    plt.legend(['D_1+','T_1+','He_1+','He_2+','Ne_1+','Ne_2+','Ne_3+','Ne_4+','Ne_5+','Ne_6+','Ne_7+','Ne_8+','Ne_9+','Ne_10+'],loc=1)
    plt.savefig('totalSputt.png')
    heTotal = np.sum(normSputtCont[:,2:4],axis=1)
    beTotal = np.sum(normSputtCont[:,4:8],axis=1)
    neTotal = np.sum(normSputtCont[:,8:18],axis=1)
    plt.close()
    plt.plot(RmRs,normSputtCont[:,0])
    plt.plot(RmRs,normSputtCont[:,1])
    plt.plot(RmRs,heTotal)
    plt.plot(RmRs,beTotal)
    plt.plot(RmRs,neTotal)
    plt.legend(['D','T','He','Be','Ne'])
    plt.savefig('totalSputt2.png')
    specArray = ['D_1+','T_1+','He_1+','He_2+','Be_1+','Be_2+','Be_3+','Be_4+','Ne_1+','Ne_2+','Ne_3+','Ne_4+','Ne_5+','Ne_6+','Ne_7+','Ne_8+','Ne_9+','Ne_10+']
    
    np.savetxt('SpylPerSpecies.dat', normSputtCont,header=" ".join(specArray))
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
    plot_hpic_ieadDavide(solps_path='/project/projectdirs/atom/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180905/solpsTarg.txt',HpicDataFolder = '/project/projectdirs/atom/atom-install-edison/GITR/data/ITER/mq4/hpicdata_solps_DT_20180905/hpicdata0005')

    #plot_hpic_ieadDavide()
    #nE,nA,nLocations,nSpecies,Z,A,dens,flux,gridE,gridA,IEAD=readHpic()
    #readFtridynSelf()
    #readFtridynBackground()
    #computeSputtYld()
    #printBackgroundDist()
    #plot_hpic_iead()
