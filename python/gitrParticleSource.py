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
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import gitr
import scipy.interpolate as scii
import netCDF4
import time
import solps

def particleSource(geomFile='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/iter/iter_milestone/3d/input/iterRefinedTest.cfg',spylFile = 'SpylPerSpecies.dat',solps_path='solpsTarg.txt',nParticles = 100000,generateFile=1):
    spyl = np.loadtxt(spylFile)
    print('spyl',spyl.shape)
    print('spyl',spyl)
    SOLPS = np.loadtxt(solps_path, dtype='float',skiprows=1,delimiter=' ')
    print('solps',SOLPS.shape)
    spFlux = 0*spyl
    rmrs = SOLPS[:,0]
    r = SOLPS[:,1]
    z = SOLPS[:,2]
    fluxInds = [28,30,32,33,35,36,37,38,40,41,42,43,44,45,46,47,48,49]
    #for i in range(len(rmrs)):
    fi  = abs(SOLPS[:,fluxInds]);
	#spFlux[i,:] = np.multiply(abs(fi),spyl[i,:])
    spFlux = np.reshape(np.sum(fi,axis=1),(36,1))*spyl
    #plt.close()
    #plt.plot(rmrs,spFlux[:,0:8])
    #plt.plot(rmrs,spFlux[:,8:16], linestyle='dashed')
    #plt.plot(rmrs,spFlux[:,16:19], linestyle='dotted')
    #plt.title('Eroded Flux By Background Species')
    #plt.ylabel('W Flux [m-2s-1]')
    #plt.xlabel('R-Rsep [m]')
    #plt.legend(['D_1+','T_1+','He_1+','He_2+','Be_1+','Be_2+','Be_3+','Be_4+','Ne_1+','Ne_2+','Ne_3+','Ne_4+','Ne_5+','Ne_6+','Ne_7+','Ne_8+','Ne_9+','Ne_10+'],loc=1)
    #plt.savefig('totalSputt.png')
    #plt.show()
    heTotal = np.sum(spFlux[:,2:4],axis=1)
    beTotal = np.sum(spFlux[:,4:8],axis=1)
    neTotal = np.sum(spFlux[:,8:18],axis=1)
    #plt.close()
    #plt.plot(rmrs,spFlux[:,0])
    #plt.plot(rmrs,spFlux[:,1])
    #plt.plot(rmrs,heTotal)
    #plt.plot(rmrs,beTotal)
    #plt.plot(rmrs,neTotal)
    #plt.title('Eroded Flux By Background Species')
    #plt.ylabel('W Flux [m-2s-1]')
    #plt.xlabel('R-Rsep [m]')
    #plt.legend(['D','T','He','Be','Ne'])
    #plt.savefig('totalSputtFlux2.png')
    filename = 'bField.nc'
    ncFile = netCDF4.Dataset(filename,"r")
    nRb = len(ncFile.dimensions['nR'])
    nZb = len(ncFile.dimensions['nZ'])
    rB = np.array(ncFile.variables['r'])
    zB = np.array(ncFile.variables['z'])
    br = np.reshape(np.array(ncFile.variables['br']),(nZb,nRb))
    bt = np.reshape(np.array(ncFile.variables['bt']),(nZb,nRb))
    bz = np.reshape(np.array(ncFile.variables['bz']),(nZb,nRb))
    
    totalFlux = spFlux[:,0]+spFlux[:,1]+heTotal+neTotal
    print(('Total flux',totalFlux))
    #print('Total flux',totalFlux.shape())
    thisRmrs = scii.interpn([z],rmrs, [-4.3])
    localFlux = scii.interpn([rmrs],totalFlux,[thisRmrs])
    print(('Local flux',localFlux))
    localFlux = scii.interpn([z],totalFlux,[-4.3])
    print(('Local flux',localFlux))
    x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray,surf,a,b,c,d,plane_norm=gitr.read3dGeom(filename=geomFile)
    rBound = 5.554-0.018
    r1 = np.sqrt(np.multiply(x1,x1) + np.multiply(y1,y1)) 
    r2 = np.sqrt(np.multiply(x2,x2) + np.multiply(y2,y2)) 
    r3 = np.sqrt(np.multiply(x3,x3) + np.multiply(y3,y3)) 
    Surf = np.where((surf>0) & (r1>=rBound) & (r2>=rBound) & (r3>=rBound) )
    xc = (x1[Surf[0]]+x2[Surf[0]]+x3[Surf[0]])/3.0
    yc = (y1[Surf[0]]+y2[Surf[0]]+y3[Surf[0]])/3.0
    zc = (z1[Surf[0]]+z2[Surf[0]]+z3[Surf[0]])/3.0
    surfaceFluxI = scii.interp1d(z,totalFlux,bounds_error=False,fill_value=0.0)
    surfaceFlux = surfaceFluxI(zc)
    surfaceParticleRate =np.multiply(surfaceFlux,area[Surf[0]])
    particleSurfs = np.where(surfaceParticleRate > 0.0)
    erodedSurfs = Surf[0][particleSurfs[0]]
    print("number of w surfaces",len(particleSurfs[0]))
    print("number of particle surfaces",len(erodedSurfs))
    particleCDF = np.cumsum(surfaceParticleRate[particleSurfs[0]])
    totalParticleRate = particleCDF[-1]
    print(('Total Particle Rate', totalParticleRate))
    particleCDF = 1/particleCDF[-1]*particleCDF
    print(('particleCDF',particleCDF))
    print(('particleCDF len',len(particleCDF)))
    #particleCDFmat = np.tile(particleCDF,(nParticles,1))
    #print('particleCDFmat',particleCDFmat.shape)
    rand1 = np.random.rand(nParticles)
    mins = np.zeros(nParticles)
    
    if generateFile:
        for i in range(nParticles):
            diff = particleCDF - rand1[i]
            diff[diff<0.0] = 100
            mins[i] = int(np.argmin(diff))
        print(('mins',mins))
        filename = 'input/angularDists.nc'
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
        #plt.plot(Egrid,Edist)
        #plt.savefig('image1.png')
        #plt.plot(thetaGrid,thetaDist)
        #plt.savefig('image2.png')
        #plt.plot(phiGrid,phiDist)
        #plt.savefig('image3.png')
        #plt.close()
        eCDF = np.cumsum(Edist,axis=0)
        eCDF = eCDF/eCDF[-1]
        phiCDF = np.cumsum(phiDist,axis=0)
        phiCDF = phiCDF/phiCDF[-1]
        thetaCDF = np.cumsum(thetaDist,axis=0)
        thetaCDF = thetaCDF/thetaCDF[-1]
        print(('ecdf sshape',eCDF.shape))
        #plt.plot(Egrid,eCDF)
        #plt.savefig('image12.png')
        #plt.plot(thetaGrid,thetaCDF)
        #plt.savefig('image22.png')
        #plt.plot(phiGrid,phiCDF)
        #plt.savefig('image32.png')
        #plt.close()
        rand2=np.random.rand(nParticles)
        rand3=np.random.rand(nParticles)
        rand4=np.random.rand(nParticles)
        #plt.hist(pE)
        #plt.savefig('image4.png')
        #plt.close()
        #plt.hist(pPhi)
        #plt.savefig('image5.png')
        #plt.close()
        #plt.hist(pTheta)
        #plt.savefig('image6.png')
        xx = np.zeros(nParticles)
        yy=np.zeros(nParticles)
        zz=np.zeros(nParticles)
        vx = np.zeros(nParticles)
        vy=np.zeros(nParticles)
        vz=np.zeros(nParticles)
        buff=1e-5
        pE = 4
        pVel = np.sqrt(2*pE*1.602e-19/184/1.66e-27)
        print(('particle velocity is ',pVel))
        pVel =np.zeros(nParticles,dtype=int)
        pPhi =np.zeros(nParticles)
        pTheta =np.zeros(nParticles)
        start = time.time()
        for i in range(nParticles):
            #print('i',i)
            #print('mins[i]',mins[i])
            #print('erodedSurfs[mins[i]]',erodedSurfs[int(mins[i])])
            minsi = int(mins[i])
            xyz = sampleTriangle(x = np.array([x1[erodedSurfs[minsi]],x2[erodedSurfs[minsi]],x3[erodedSurfs[minsi]]]),
            y = np.array([y1[erodedSurfs[minsi]],y2[erodedSurfs[minsi]],y3[erodedSurfs[minsi]]]),
            z = np.array([z1[erodedSurfs[minsi]],z2[erodedSurfs[minsi]],z3[erodedSurfs[minsi]]]))
            xx[i]=(xyz[0] - a[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*buff)
            yy[i]=(xyz[1] - b[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*buff)
            zz[i]=(xyz[2] - c[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*buff)
        
            rzDiff = z-zz[i]
            rmrsInd = int(np.argmin(abs(rzDiff)))
            #pEint = scii.interp1d(eCDF[:,rmrsInd],Egrid,bounds_error=False,fill_value=4.0)
            eCDFdiff = eCDF[:,rmrsInd] - rand2[i]
            eCDFInd = int(np.argmin(abs(eCDFdiff)))
            pE = Egrid[eCDFInd]
            #pE = pEint(rand2[i])
            if pE > 20.0:
                pE = 6.0
            #pE = 250.0
            pVel[i] = np.sqrt(2*pE*1.602e-19/184/1.66e-27)
            #pPhiint = scii.interp1d(phiCDF[:,rmrsInd],phiGrid,bounds_error=False,fill_value=0.0)
            #pPhi = pPhiint(rand3[i])
            phiCDFdiff = phiCDF[:,rmrsInd] - rand3[i]
            phiCDFInd = int(np.argmin(abs(phiCDFdiff)))
            pPhi[i] = phiGrid[phiCDFInd]

            #pThetaint = scii.interp1d(thetaCDF[:,rmrsInd],thetaGrid,bounds_error=False,fill_value=90.0)
            #pTheta = pThetaint(rand4[i])
            thetaCDFdiff = thetaCDF[:,rmrsInd] - rand4[i]
            thetaCDFInd = int(np.argmin(abs(thetaCDFdiff)))
            pTheta[i] = thetaGrid[thetaCDFInd]
            #pTheta = 90
            #pPhi=0.0
        vxSampled = pVel*np.sin(pPhi*3.1415/180)*np.cos(3.1415*pTheta/180.0);
        vySampled = pVel*np.sin(pPhi*3.1415/180)*np.sin(3.1415*pTheta/180.0);
        vzSampled = pVel*np.cos(pPhi*3.1415/180);
        thisR = np.sqrt(xx*xx + yy*yy)
        rzArray = np.zeros((nParticles,2))
        rzArray[:,1] = thisR
        rzArray[:,0] = zz
        locBr = scii.interpn([zB,rB],br, rzArray)
        locBt = scii.interpn([zB,rB],bt, rzArray)
        locBz = scii.interpn([zB,rB],bz, rzArray)
        theta = np.arctan2(yy,xx);   
        localfield0 = np.cos(theta)*locBr - np.sin(theta)*locBt;
        localfield1 = np.sin(theta)*locBr + np.cos(theta)*locBt;
        for i in range(nParticles):
            minsi = int(mins[i])
            surfNormx = ( - a[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]])
            surfNormy = ( - b[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]])
            surfNormz = ( - c[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]])
            bbb = [surfNormx, surfNormy, surfNormz]
            #bbb = [-surfNormx, -surfNormy, -surfNormz]
            aaa = [localfield0[i], localfield1[i], locBz[i]]
            #bbb = [0, 0, 1]
            #aaa = [0, 1, 0]
            #print('aaa',aaa)
            #print('bbb',bbb)
            
            #surfParallelX = x2[erodedSurfs[minsi]]-x1[erodedSurfs[minsi]] 
            #surfParallelY = y2[erodedSurfs[minsi]]-y1[erodedSurfs[minsi]] 
            #surfParallelZ = z2[erodedSurfs[minsi]]-z1[erodedSurfs[minsi]] 
            #surfParallelNorm = sqrt(surfParallelX*surfParallelX + surfParallelY*surfParallelY + surfParallelZ*surfParallelZ)
            #surfParallelX = surfParallelX/surfParallelNorm
            #surfParallelY = surfParallelY/surfParallelNorm
            #surfParallelZ = surfParallelZ/surfParallelNorm
            surfPar = np.cross(bbb,np.cross(aaa,bbb))
            surfParNorm = np.linalg.norm(surfPar)
            surfPar = surfPar/surfParNorm
            surfParX = np.cross(surfPar,bbb)
            #print('surfParx',surfParX)
            #print('surfPar',surfPar)
            vx[i]=(surfParX[0]*vxSampled[i] + surfPar[0]*vySampled[i] + bbb[0]*vzSampled[i] )
            vy[i]=(surfParX[1]*vxSampled[i] + surfPar[1]*vySampled[i] + bbb[1]*vzSampled[i])
            vz[i]=(surfParX[2]*vxSampled[i] + surfPar[2]*vySampled[i] + bbb[2]*vzSampled[i])
            #vx.append( - a[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*pVel)
            #vy.append( - b[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*pVel)
            #vz.append( - c[erodedSurfs[minsi]]/plane_norm[erodedSurfs[minsi]]*pVel)
            
        #xx = x1[erodedSurfs[mins]]
        #yy = y1[erodedSurfs[mins]]
        #zz = z1[erodedSurfs[mins]]
        #print('xx',xx)
        #print('yy',yy)
        #print('zz',zz)
        end = time.time()
        print(('Time in loop',end - start))
        x = []
        y = []
        z = []
        for j in range(len(erodedSurfs)):
            i = erodedSurfs[j]
            x.append(x1[i])
            x.append(x2[i])
            x.append(x3[i])
            y.append(y1[i])
            y.append(y2[i])
            y.append(y3[i])
            z.append(z1[i])
            z.append(z2[i])
            z.append(z3[i])
        x=np.array(x)
        y=np.array(y)
        z=np.array(z)
        print('mininum z',z.min())
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.scatter(xc, yc, zc, c=surfaceParticleRate, marker='o')
        #ax.scatter(xx, yy, zz, c='r', marker='o')
        ##ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
        #plt.show()
        #plt.savefig('surfaces.png')
        #plt.close()
        rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
        npp = rootgrp.createDimension("nP", nParticles)
        xxx = rootgrp.createVariable("x","f8",("nP"))
        yyy = rootgrp.createVariable("y","f8",("nP"))
        zzz = rootgrp.createVariable("z","f8",("nP"))
        vxx = rootgrp.createVariable("vx","f8",("nP"))
        vyy = rootgrp.createVariable("vy","f8",("nP"))
        vzz = rootgrp.createVariable("vz","f8",("nP"))
        xxx[:] = xx
        yyy[:] = yy
        zzz[:] = zz
        vxx[:] = vx
        vyy[:] = vy
        vzz[:] = vz
        rootgrp.close()
    return totalParticleRate 
def sampleTriangle(x=np.array([1, 1.2, 1.3]),y=np.array([0, 0.1, -0.1]),z=np.array([1.1, 1.05, 0.92])):
    x_transform = x - x[0]
    y_transform = y - y[0]
    z_transform = z - z[0]
    v1 = np.array([x_transform[1], y_transform[1], z_transform[1]]);
    v2 = np.array([x_transform[2], y_transform[2], z_transform[2]]);
    v12 = v2 - v1;
    normalVec = np.cross(v1,v2);

    a1 = np.random.rand(1);
    a2 = np.random.rand(1);

    samples = a1*v1 + a2*v2;
    samples2 = samples - v2
    samples12 = samples - v1
    v1Cross = np.cross(v1,samples)
    v2=-v2
    v2Cross = np.cross(v2,samples2)
    v2 = -v2;
    v12Cross = np.cross(v12,samples12)
    
    v1CD = np.dot(normalVec,v1Cross)
    v2CD = np.dot(normalVec,v2Cross)
    v12CD = np.dot(normalVec,v12Cross)
    
    inside = np.absolute(np.sign(v1CD) + np.sign(v2CD) + np.sign(v12CD));
    isInside = 0
    if inside ==3:
        isInside=1
    else:
        dAlongV1 = np.dot(v1,samples)
        dAlongV2 = np.dot(v2,samples)
    
        dV1 = np.linalg.norm(v1);
        dV2 = np.linalg.norm(v2);
        halfdV1 = 0.5*dV1;
        halfdV2 = 0.5*dV2;
    
        samples = [-(samples[0] - 0.5*v1[0])+0.5*v1[0],
                   -(samples[1] - 0.5*v1[1])+0.5*v1[1],
                   -(samples[2] - 0.5*v1[2])+0.5*v1[2]];
        samples = [(samples[0] + v2[0]),
	           (samples[1] + v2[1]),
	           (samples[2] + v2[2])];
		 
    samples[0] = samples[0]+ x[0];
    samples[1] = samples[1]+ y[0];
    samples[2] = samples[2]+ z[0];
    return samples
def particleSource2d(nParticles = int(1e3),spylFile = 'Yields.txt', coordsFile='right_target_coordinates.txt',edist_file = 'Edists.txt',adist_file = 'Adists.txt'):
    y_f = np.loadtxt(spylFile, dtype='float',skiprows=1,delimiter=' ')
    eDists = np.loadtxt(edist_file, dtype='float',skiprows=1,delimiter=' ')
    aDists = np.loadtxt(adist_file, dtype='float',skiprows=1,delimiter=' ')
    coords = np.loadtxt(coordsFile, dtype='float',skiprows=1,delimiter=' ')
    x1, x2, z1, z2, length, Z, slope, inDir = gitr.plot2dGeom('gitr_geometry.cfg')
    print('shapes',y_f.shape,coords.shape,eDists.shape)

    #calcualte erosion flux
    yf_shape = y_f.shape
    nSpec = int(0.5*yf_shape[1])

    sputt_flux = np.zeros((yf_shape[0],nSpec+1))

    for i in range(nSpec):
        sputt_flux[:,i] = np.multiply(y_f[:,i],y_f[:,i+nSpec])
        sputt_flux[:,-1] = sputt_flux[:,-1] + sputt_flux[:,i]

    #calculate ribbon area
    gitr_inds = [int(i) for i in coords[0:-1,3]]
    print('gitr_inds',gitr_inds)
    r1 = x1[gitr_inds] #coords[:,4]
    r2 = x2[gitr_inds] #coords[:,5]
    z1 = z1[gitr_inds] #coords[:,6]
    z2 = z2[gitr_inds] #coords[:,7]
    slope = slope[gitr_inds] #coords[:,7]
    inDir = inDir[gitr_inds] #coords[:,7]
    length = length[gitr_inds] #coords[:,7]
    print('slope',slope)
    
    area = np.pi*(r1+r2)*np.sqrt(np.power(r1-r2,2) + np.power(z1 - z2,2))
    print('area',area)
    
    print('len area',len(area[0:-1]))
    print('len sputt flux', len(sputt_flux[:,-1]))    

    particles_per_second = np.multiply(np.abs(sputt_flux[:,-1]),area)

    print('sputt flux', sputt_flux)    
    print('area',area)
    print('pps',particles_per_second)
    pps_cdf = np.cumsum(particles_per_second)
    pps_sum = pps_cdf[-1]
    pps_cdf = pps_cdf/pps_cdf[-1]
    pps_cdf = np.transpose(np.matlib.repmat(pps_cdf,nParticles,1))
    rand1 = np.random.rand(nParticles)
    rand1 = np.matlib.repmat(rand1,yf_shape[0],1)

    print('rand1 shape',rand1.shape)

    diff = pps_cdf - rand1
    diff[diff<0.0] = 100.0
    mins = np.argmin(diff,axis=0)
    
    print('len mins',len(mins))
    print('mins',mins)
    print('diff',diff[mins,range(nParticles)])
    print('scale',pps_cdf[-1]/particles_per_second[mins])

    r_sample = r1[mins] + (r2[mins] - r1[mins])*diff[mins,range(nParticles)]/particles_per_second[mins]*pps_sum
    z_sample = z1[mins] + (z2[mins] - z1[mins])*diff[mins,range(nParticles)]/particles_per_second[mins]*pps_sum
    print('rsamp',r_sample)
    print('zsamp',z_sample)
    plt.close()
    plt.plot(r1,z1)
    plt.scatter(r_sample,z_sample)
    plt.title('Outer divertor')
    plt.ylabel('z [m]')
    plt.xlabel('r [m]')
    plt.savefig('particleScatter.png')
    rPar = np.divide((r2 - r1),length)
    zPar = np.divide((z2 - z1),length)
    #print('should be ones', np.sqrt(np.multiply(rPar,rPar) + np.multiply(zPar,zPar)))
    parVec = np.zeros([len(gitr_inds),3])
    parVec[:,0] = rPar
    parVec[:,2] = zPar
    rPerp = 0*slope;
    zPerp = 0*slope;
    for i in range(0,len(slope)):
         if slope[i]==0:
             perpSlope = 1.0e12
         else:
             perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

         rPerp[i] = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
         zPerp[i] = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp[i]*rPerp[i]);
    #perpSlope = -np.sign(slope)/np.abs(slope);
    #rPerp = -inDir/np.sqrt(perpSlope*perpSlope+1);
    #zPerp = -inDir*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
    perpVec = np.zeros([len(gitr_inds),3])
    perpVec[:,0] = rPerp
    perpVec[:,2] = zPerp
    print('perpVec',perpVec)
    
    surface_buffer = 1.0e-5
    r_sample = r_sample + surface_buffer*rPerp[mins]
    z_sample = z_sample + surface_buffer*zPerp[mins]

    particle_energy = np.zeros(nParticles)
    particle_angle = np.zeros(nParticles)
    mins_array = np.array(mins)
    plt.close()
    print('len gitr inds',len(gitr_inds))
    for i in range(len(gitr_inds)):
        e_cdf = np.cumsum(eDists[i,:])
        plt.plot(eDists[i,:])
        e_sum = e_cdf[-1]
        if (e_sum > 0.0):
            e_cdf = e_cdf/e_sum
            if (e_cdf[0] != 0.0):
                e_cdf[0] = 0.0
            where_mins = np.where(mins_array == i)
            print('where mins',where_mins)
            if (len(where_mins[0]) > 0):
                these_energies = gitr.interp_1d(e_cdf,np.linspace(0.0,40.0,40),np.random.rand(len(where_mins[0])),default_value = 5.0)
                particle_energy[where_mins[0]] = these_energies
    plt.savefig('edist.png')
    plt.close()
    plt.hist(particle_energy)
    plt.savefig('pE.png')
    plt.close()
    
    particle_angle = np.zeros(nParticles)
    mins_array = np.array(mins)
    plt.close()
    for i in range(len(gitr_inds)):
        a_cdf = np.cumsum(aDists[i,:])
        plt.plot(aDists[i,:])
        a_sum = a_cdf[-1]
        if (a_sum > 0.0):
            a_cdf = a_cdf/a_sum
            if (a_cdf[0] != 0.0):
                a_cdf[0] = 0.0
            where_mins = np.where(mins_array == i)
            print('where mins',where_mins)
            if (len(where_mins[0]) > 0):
                these_angles = gitr.interp_1d(a_cdf,np.linspace(0.0,89.9,50),np.random.rand(len(where_mins[0])),default_value = 45.0)
                particle_angle[where_mins[0]] = these_angles
    plt.savefig('adist.png')
    plt.close()
    plt.hist(particle_angle)
    plt.savefig('pA.png')
    plt.close()

    randDistvTheta = np.random.rand(nParticles)*2*np.pi
    vr = np.sqrt(2*particle_energy*1.602e-19/184/1.66e-27);
    #vx = np.multiply(vr,np.multiply(np.cos(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vy = np.multiply(vr,np.multiply(np.sin(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vz = np.multiply(vr,np.cos(np.pi/180.0*particle_angle));
    vx = 0*vr;
    vy = 0*vr;
    vz = vr;
    parVec2 = np.cross(parVec,perpVec)
    vxP = 0*vr;
    vyP = 0*vr;
    vzP = 0*vr;

    
    for i in range(nParticles):
        vxP[i] = parVec[mins[i],0]*vx[i] + parVec2[mins[i],0]*vy[i] + perpVec[mins[i],0]*vz[i];
        vyP[i] = parVec[mins[i],1]*vx[i] + parVec2[mins[i],1]*vy[i] + perpVec[mins[i],1]*vz[i];
        vzP[i] = parVec[mins[i],2]*vx[i] + parVec2[mins[i],2]*vy[i] + perpVec[mins[i],2]*vz[i];
        if z_sample[i] < -3.8 :
            print('r,z',r_sample[i],z_sample[i])
            print('perpVec xyz',perpVec[mins[i],0],perpVec[mins[i],1],perpVec[mins[i],2])
            print('v xyz',vxP[i],vyP[i],vzP[i])

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nParticles)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = r_sample
    yyy[:] = 0*r_sample
    zzz[:] = z_sample
    vxx[:] = vxP
    vyy[:] = vyP
    vzz[:] = vzP
    rootgrp.close()

def particleSource2d_west(nParticles = int(1e3), \
                            geom='/Users/Alyssa/Dev/GITR/west/helium/input/gitrGeometry.cfg', \
                            targFile='/Users/Alyssa/Dev/solps-iter-data/build/rightTargOutput', \
                            coordsFile='/Users/Alyssa/Dev/GITR/west/helium/output/right_target_coordinates.txt'):

    x1, x2, z1, z2, length, Z, slope, inDir = gitr.plot2dGeom(geom)
    coords = np.loadtxt(coordsFile, dtype='float',skiprows=1,delimiter=' ')
    gitr_inds = [int(i) for i in coords[0:-1,3]]
    print('gitr_inds',gitr_inds)
    r1 = x1[gitr_inds] #coords[:,4]
    r2 = x2[gitr_inds] #coords[:,5]
    z1 = z1[gitr_inds] #coords[:,6]
    z2 = z2[gitr_inds] #coords[:,7]
    slope = slope[gitr_inds] #coords[:,7]
    inDir = inDir[gitr_inds] #coords[:,7]
    length = length[gitr_inds] #coords[:,7]
    area = np.pi*(r1+r2)*np.sqrt(np.power(r1-r2,2) + np.power(z1 - z2,2))

    r, z, ti, ni, flux, te, ne = solps.read_target_file(targFile)

    ion_flux = flux[:,1:]
    yf = np.hstack((0.1+0*ion_flux,ion_flux))
    nSpec = int(0.5*yf.shape[1])
    print('nSpec',nSpec)

    sputt_flux = np.zeros((yf.shape[0],nSpec+1))

    for i in range(nSpec):
        sputt_flux[:,i] = np.multiply(yf[:,i],yf[:,i+nSpec])
        sputt_flux[:,-1] = sputt_flux[:,-1] + sputt_flux[:,i]

    buff = np.array([[0,0,0]])
    sputt_flux = np.vstack((buff,sputt_flux,buff))
    print('len sputt flux', len(sputt_flux[:,-1]))

    particles_per_second = np.multiply(np.abs(sputt_flux[:,-1]),area)

    print('sputt flux', sputt_flux)
    print('area',area)
    print('pps',particles_per_second)
    pps_cdf = np.cumsum(particles_per_second)
    pps_sum = pps_cdf[-1]
    pps_cdf = pps_cdf/pps_cdf[-1]
    pps_cdf = np.transpose(np.matlib.repmat(pps_cdf,nParticles,1))
    rand1 = np.random.rand(nParticles)
    rand1 = np.matlib.repmat(rand1,sputt_flux.shape[0],1)

    print('rand1 shape',rand1.shape)

    diff = pps_cdf - rand1
    diff[diff<0.0] = 100.0
    mins = np.argmin(diff,axis=0)

    print('len mins',len(mins))
    print('mins',mins)
    print('diff',diff[mins,range(nParticles)])
    print('scale',pps_cdf[-1]/particles_per_second[mins])

    r_sample = r1[mins] + (r2[mins] - r1[mins])*diff[mins,range(nParticles)]/particles_per_second[mins]*pps_sum
    z_sample = z1[mins] + (z2[mins] - z1[mins])*diff[mins,range(nParticles)]/particles_per_second[mins]*pps_sum
    #print('rsamp',r_sample)
    #print('zsamp',z_sample)
    plt.close()
    plt.plot(r1,z1)
    plt.scatter(r_sample,z_sample,alpha=0.1)
    plt.title('Outer divertor')
    plt.ylabel('z [m]')
    plt.xlabel('r [m]')
    plt.savefig('particleScatter.png')

    plt.close()
    plt.hist(r_sample)
    plt.savefig('hist.png')

    rPar = np.divide((r2 - r1),length)
    zPar = np.divide((z2 - z1),length)
    #print('should be ones', np.sqrt(np.multiply(rPar,rPar) + np.multiply(zPar,zPar)))
    parVec = np.zeros([len(gitr_inds),3])
    parVec[:,0] = rPar
    parVec[:,2] = zPar

    print('inDir', inDir)
    rPerp = 0*slope;
    zPerp = 0*slope;
    for i in range(0,len(slope)):
         if slope[i]==0:
             perpSlope = 1.0e12
         else:
             perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

         rPerp[i] = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
         zPerp[i] = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp[i]*rPerp[i]);
    #perpSlope = -np.sign(slope)/np.abs(slope);
    #rPerp = -inDir/np.sqrt(perpSlope*perpSlope+1);
    #zPerp = -inDir*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
    perpVec = np.zeros([len(gitr_inds),3])
    perpVec[:,0] = rPerp
    perpVec[:,2] = zPerp
    print('perpVec',perpVec)

    #moves particles 10um away from surface
    surface_buffer = 1.0e-5
    r_sample = r_sample + surface_buffer*rPerp[mins]
    z_sample = z_sample + surface_buffer*zPerp[mins]

    particle_energy = 4*np.ones(nParticles)
    particle_angle = np.zeros(nParticles)

    randDistvTheta = np.random.rand(nParticles)*2*np.pi
    vr = np.sqrt(2*particle_energy*1.602e-19/184/1.66e-27);
    #vx = np.multiply(vr,np.multiply(np.cos(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vy = np.multiply(vr,np.multiply(np.sin(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vz = np.multiply(vr,np.cos(np.pi/180.0*particle_angle));
    vx = 0*vr;
    vy = 0*vr;
    vz = vr;
    parVec2 = np.cross(parVec,perpVec)
    vxP = 0*vr;
    vyP = 0*vr;
    vzP = 0*vr;

    for i in range(nParticles):
        vxP[i] = parVec[mins[i],0]*vx[i] + parVec2[mins[i],0]*vy[i] + perpVec[mins[i],0]*vz[i];
        vyP[i] = parVec[mins[i],1]*vx[i] + parVec2[mins[i],1]*vy[i] + perpVec[mins[i],1]*vz[i];
        vzP[i] = parVec[mins[i],2]*vx[i] + parVec2[mins[i],2]*vy[i] + perpVec[mins[i],2]*vz[i];
        #if z_sample[i] < -3.8 :
            #print('r,z',r_sample[i],z_sample[i])
            #print('perpVec xyz',perpVec[mins[i],0],perpVec[mins[i],1],perpVec[mins[i],2])
            #print('v xyz',vxP[i],vyP[i],vzP[i])

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nParticles)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = r_sample
    yyy[:] = 0*r_sample
    zzz[:] = z_sample
    vxx[:] = vxP
    vyy[:] = vyP
    vzz[:] = vzP
    rootgrp.close()

if __name__ == "__main__":
    particleSource2d_west()
    #x = []
    #y=[]
    #z=[]
    #for i in range(100):
    #    xyz = sampleTriangle()
    #    x.append(xyz[0])
    #    y.append(xyz[1])
    #    z.append(xyz[2])
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax.scatter(x, y, z, c='r', marker='o')
    ##ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
    #plt.show()
