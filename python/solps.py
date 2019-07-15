import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as ml
import gitr
import scipy.interpolate as scii
import netCDF4
def readEquilibrium(filename='/global/homes/t/tyounkin/atomIPS/atom-install-edison/solps-iter-data/Baseline2008-li0.70.x4.equ',geometryFile='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/iter/iter_milestone/2d/input/iter2dRefinedOuterTarget.cfg'):
    print('Reading B Equilibrium file %s and making use of GITR geometry file %s' %(filename,geometryFile))
    rr=0
    zz=0
    pp=0
    r = []
    z = []
    psi = []
    with open(filename) as openfileobject:
        for line in openfileobject:
                #print line
            if line:
    	    #    print 'empty line'
    	    #else:
                    l = line.split()
                    if not l:
                        k=1
    		    #print 'empty line'
                    else:
    		    #print l
                        if (l[0] == "jm" and l[1]=="="):
                            jm = int(l[2])
                        if (l[0] == "km" and l[1]=="="):
                            km = int(l[2])
                        if (l[0] == "psib" and l[1]=="="):
                            psib = float(l[2])
                        if (l[0] == "btf" and l[1]=="="):
                            btf = float(l[2])
                        if (l[0] == "rtf" and l[1]=="="):
                            rtf = float(l[2])
                        if (l[0] == "((psi(j,k)-psib,j=1,jm),k=1,km)" or pp==1):
                            if pp:
                                ll=[float(i) for i in l]
                                psi = np.concatenate([psi, ll])
                            zz=0;
                            pp=1;
                        if (l[0] == "z(1:km);" or zz==1):
                            if zz:
                                ll=[float(i) for i in l]
                                z = np.concatenate([z, ll])
                            rr=0;
                            zz=1;
                        if (l[0] == "r(1:jm);" or rr==1):
                            if rr:
                                ll=[float(i) for i in l]
                                r = np.concatenate([r, ll])
                            rr=1;
    
    
    print ('Equ data dimensions %i by %i ' %(jm,km)	)
    psi = np.reshape(psi,[len(z),len(r)])
    #print(psi.shape)
    plt.pcolor(r, z, psi)
    plt.title('pcolor')
    # set the limits of the plot to the limits of the data
    plt.axis([r.min(), r.max(), z.min(), z.max()])
    plt.colorbar()
    print( 'Saving psi function as psi.png ')
    plt.savefig('psi.png')
    plt.close()
    plt.contour(r,z,psi,100)
    plt.title('contour')
    # set the limits of the plot to the limits of the data
    plt.axis([r.min(), r.max(), z.min(), z.max()])
    plt.colorbar()
    gitr.plot2dGeom(geometryFile)
    print ('Saving psi contour as psiContour.png ')
    plt.savefig('psiContour.png')
    print('About to take gradient')
    print(r.shape)
    print(z.shape)
    print(psi.shape)
    [gradz,gradr] = np.gradient(np.array(psi),z[1]-z[0],r[1]-r[0]) #z,r) sometimes this doesn't work
    print(gradr.shape)
    print(gradz.shape)
    br = -gradz/r
    print(br.shape)
    bz =gradr/r
    print(bz.shape)
    plt.close()
    print(r)
    print(z)
    plt.pcolor(r,z,br)
    plt.colorbar()
    print ('Saving br profile as br.png ')
    plt.savefig('br.png')
    plt.close()
    plt.pcolor(r,z,bz)
    plt.colorbar()
    print ('Saving br profile as br.png ')
    plt.savefig('bz.png')
    Bp = np.sqrt(np.multiply(br,br) + np.multiply(bz,bz))
    bt = ml.repmat(btf*rtf/r,len(z),1)
    bt = np.reshape(bt,[len(z),len(r)])
    plt.close()
    plt.pcolor(r,z,bt)
    plt.colorbar()
    print( 'Saving br profile as br.png ')
    plt.savefig('bt.png')
    plt.close()
    #print 'plotting psi'
    #plt.contour(r,z,psi,np.linspace(-2,2,20),linewidths=0.3)
    ##plt.axis.set_aspect('equal')
    ##plt.axis.autoscale(tight=True)
    ##plt.axis({'equal','tight'})
    ##plt.axis('equal')
    ##plt.axis('tight')
    ##plt.autoscale(enable=True, axis='x', tight=True)
    #plt.colorbar()
    #gitr.plot2dGeom(geometryFile)
    #plt.axis('equal')
    #plt.axis('tight')
    #plt.savefig('psiContour.png')
    rootgrp = netCDF4.Dataset("bField.nc", "w", format="NETCDF4")
    nrr = rootgrp.createDimension("nR", len(r))
    nzz = rootgrp.createDimension("nZ", len(z))
    brr = rootgrp.createVariable("br","f8",("nZ","nR"))
    btt = rootgrp.createVariable("bt","f8",("nZ","nR"))
    bzz = rootgrp.createVariable("bz","f8",("nZ","nR"))
    rr = rootgrp.createVariable("r","f8",("nR"))
    zz = rootgrp.createVariable("z","f8",("nZ"))
    brr[:] = br
    btt[:] = bt
    bzz[:] = bz
    rr[:] = r
    zz[:] = z
    rootgrp.close()
    print( 'finishing read equilibirium')
    return r,z,br, bz, bt,psi
def findStrikepoint(x1,x2,z1,z2,length,r,z,psi,rmin=5.55,rmax=6.226,zmin=-4.6,zmax=-3.238):
    outerTargetInd = np.where((x1 > rmin) & (x1 < rmax) & (x2 > rmin) & (x2 < rmax) & (z1 > zmin) & (z1 < zmax) & (z2 > zmin) & (z2 < zmax))
    #print 'outerTargetInd', outerTargetInd  
    targInds = outerTargetInd[0]
    #print 'targInds', targInds
    psiTarg = np.zeros((len(targInds),2))
    for i in range(len(targInds)):
        targNum = targInds[i]
	    #print 'loop', i, targNum

    psiTarg[i,0] = scii.interpn([r,z],psi.T, [x1[targNum],z1[targNum]])
    psiTarg[i,1] = scii.interpn([r,z],psi.T, [x2[targNum],z2[targNum]])
    sepBoundary = np.where((np.sign(psiTarg[:,0]) < 0.0) & (np.sign(psiTarg[:,1]) > 0.0))
    sepBoundary2 = targInds[sepBoundary];
    rsep = x1[sepBoundary2] + abs(psiTarg[sepBoundary,0])/(psiTarg[sepBoundary,1]-psiTarg[sepBoundary,0])*(x2[sepBoundary2]-x1[sepBoundary2])
    zsep = z1[sepBoundary2] + abs(psiTarg[sepBoundary,0])/(psiTarg[sepBoundary,1]-psiTarg[sepBoundary,0])*(z2[sepBoundary2]-z1[sepBoundary2])
    print( 'Rsep Zsep of strike point: ',rsep,zsep)
    return rsep, zsep, targInds, sepBoundary2
def interpolateBfield(r,z,br, bz, bt,psi,rTarget,zTarget,geometryFile='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/iter/iter_milestone/2d/input/iter2dRefinedOuterTarget.cfg',rmin=5.55,rmax=6.226,zmin=-4.6,zmax=-3.238):
    print ('Beginning interpolate B-field')
    x1,x2,z1,z2,length,Z = gitr.plot2dGeom(geometryFile) 
    print( 'Finding outer target strike point')
    rsep,zsep,targInds,sepBoundary = findStrikepoint(x1,x2,z1,z2,length,r,z,psi,rmin,rmax,zmin,zmax)
    nS = len(rTarget)
    #for i in range(nS)
    plt.close()
    rMrs = np.zeros((len(targInds)+1,))
    rTarg = np.zeros((len(targInds)+1,))
    zTarg = np.zeros((len(targInds)+1,))
    slope = np.zeros((len(targInds),))
    surfaceNormal = np.zeros((len(targInds),3))
    #print 'rMrs shape and rTarg.shape',rMrs.shape,rTarg.shape
    print( 'Computing R-Rsep and Bfield angle at target locations')
    rTarg[0] = x2[targInds[0]]
    zTarg[0] = z2[targInds[0]]
    for j in targInds:
        ref = j-targInds[0]+1
        rMrs[ref] = rMrs[ref-1]+length[j]
        rTarg[ref] = x2[j+1]
        zTarg[ref] = z2[j+1]

        plt.plot([x1[j],x2[j]],[z1[j],z2[j]])
    plt.text(x1[j],z1[j],str(j),size=6)
    rise = z1[j] - z2[j]
    run = x1[j] - x2[j]
    slope[ref-1] = rise/run
    intercept=0
    if(( zsep > z2[j]) and (zsep < z1[j])):
        intercept = z1[j] - x1[j]*slope[ref-1]
        xsep = (zsep-intercept)/slope[ref-1]
        dx2 = np.sqrt((xsep-x2[j])**2 + (zsep - z2[j])**2)
        dx1 = np.sqrt((xsep-x1[j])**2 + (zsep - z1[j])**2)
        rMrsShift = rMrs[ref-1] + dx2
        if(np.isinf(slope[ref-1]) or np.isnan(slope[ref-1])):
            surfaceNormal[ref-1,0] = 1
        elif (slope[ref-1] == 0):
            surfaceNormal[ref-1,2] = 1
        else:
            surfaceNormal[ref-1,0] = 1/np.sqrt(1+1/slope[ref-1]**2)
            surfaceNormal[ref-1,2] = -1/slope[ref-1]/np.sqrt(1+1/slope[ref-1]**2)
            tot = np.sqrt(surfaceNormal[ref-1,0]**2+surfaceNormal[ref-1,2]**2)
            surfaceNormal[ref-1,0]=surfaceNormal[ref-1,0]/tot
        surfaceNormal[ref-1,2]=surfaceNormal[ref-1,2]/tot
    rMrs = rMrs - float(rMrsShift)
    #print 'rMrs shape and rTarg.shape',rMrs.shape,rTarg.shape
    #print rMrs
    #print 'rMrs shape and rTarg.shape',rMrs.shape,rTarg.shape
    plt.savefig('target.png')	
    rloc = scii.interp1d(rMrs,rTarg) 
    zloc = scii.interp1d(rMrs,zTarg) 
    #print 'rloc zloc', rloc(0), zloc(0)
    startingIndex=targInds[0]
    endIndex = targInds[-1]
    surfNum = np.abs(rMrs-0).argmin() + targInds[0]
    #print 'surfNum ', surfNum, z1[surfNum],z2[surfNum]
    #print 'normal ', surfaceNormal[surfNum-startingIndex,:]
    #print 'normal ', surfaceNormal
    
    #print 'process rsep'
    bAngle = np.zeros(nS)
    bMag = np.zeros(nS)
    rSep = np.zeros(nS)
    x1Targ=x1[int(startingIndex):int(endIndex)]
    x2Targ=x2[int(startingIndex):int(endIndex)]
    z1Targ=z1[int(startingIndex):int(endIndex)]
    z2Targ=z2[int(startingIndex):int(endIndex)]
    #print 'x1Targ size shape',x1Targ.size,x1Targ.shape
    #print 'start end',startingIndex,endIndex
    for i in range(nS):
        br0 = scii.interpn([z,r],br, [zTarget[i],rTarget[i]]) 
        bz0 = scii.interpn([z,r],bz, [zTarget[i],rTarget[i]]) 
        bt0 = scii.interpn([z,r],bt, [zTarget[i],rTarget[i]])
        bMag[i] = np.sqrt(br0*br0 + bz0*bz0 + bt0*bt0)
        rMrSep,surfNum =getRsepFromRZ(x1Targ,x2Targ,z1Targ,z2Targ,slope,rMrs,rTarget[i],zTarget[i])
        rSep[i] = rMrSep
        surfNum=surfNum + startingIndex
        #print 'surfNum',surfNum
        normal = surfaceNormal[surfNum-startingIndex,:]
        bAngle[i] = 180/np.pi*np.arccos(np.dot(normal,[br0,bt0,bz0])/bMag[i])
        #print 'brtz',br0,bt0,bz0
    return rSep,bAngle,bMag
def getRsepFromRZ(x1,x2,z1,z2,slope,rMrs,r,z):
    #print 'get rsep from rz',r,z
    nS = len(x1)
    #print 'nS',nS
    surfNumber = -1
    for j in range(nS):
        if(( z >= z2[j]) and (z <= z1[j])):
            #intercept = z1[j] - x1[j]*slope[j]
            #print 'gt z2 lt z1',z2[j],z1[j],j
            #rsep = (zsep-intercept)/slope[ref-1]
            dx2 = np.sqrt((r-x2[j])**2 + (z - z2[j])**2)
            surfNumber = j
    if surfNumber ==-1:
        #print 'position out of target rz',r,z
	    rMrSep=0
    else:
        rMrSep = rMrs[surfNumber]+dx2
    return rMrSep,surfNumber
def getBfield(rTarg,zTarg,filename='/global/homes/t/tyounkin/atomIPS/atom-install-edison/solps-iter-data/Baseline2008-li0.70.x4.equ',geometryFile='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/iter/iter_milestone/2d/input/iter2dRefinedOuterTarget.cfg',rmin=5.55,rmax=6.226,zmin=-4.6,zmax=-3.238):    
    r,z,br, bz, bt,psi = readEquilibrium(filename,geometryFile)
    rSep,bAngle, bMag = interpolateBfield(r,z,br, bz, bt,psi,rTarg,zTarg,geometryFile,rmin,rmax,zmin,zmax)
    print( 'plotting bangle at target')
    plt.close()
    plt.plot(rSep,bAngle)
    plt.savefig('bangle.png')
    return rSep,bAngle,bMag
if __name__ == "__main__":   
    rTarg = np.linspace(5,6.5,100)
    zTarg=np.linspace(0,1,100)
    getBfield(rTarg,zTarg,"/Users/tyounkin/Dissertation/ITER/mq3/final/Baseline2008-li0.70.x4.equ","/Users/tyounkin/Code/gitr2/iter/iter_milestone/2d/input/iterGeom2DdirBe0.cfg")
