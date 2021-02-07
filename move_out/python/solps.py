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
from scipy.interpolate import griddata
import netCDF4
import math
import os
import io
import libconf
import solps

def readEquilibrium(filename='/Users/Alyssa/Dev/WEST/baserun/west_54034_10p2s_mag.X4.equ', \
    solps_mesh_extra = '/Users/Alyssa/Dev/WEST/baserun/mesh.extra', \
    solps_geom = '/Users/Alyssa/Dev/WEST/baserun/b2fgmtry'):

    rr=0
    zz=0
    pp=0
    r = []
    z = []
    psi = []
    with open(filename) as openfileobject:
        for line in openfileobject:
            if line:
                    l = line.split()
                    if not l:
                        k=1
                    else:
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


    #get geometry to plot on top of magnetic fields / fluxes
    solps_mesh = np.loadtxt(solps_mesh_extra)
    r_geom = solps_mesh[:, [0,2]].transpose()
    z_geom = solps_mesh[:, [1,3]].transpose()
    r_left_target,z_left_target,r_right_target,z_right_target = solps.get_target_coordinates(solps_geom)

    print ('Equ data dimensions %i by %i ' %(jm,km)	)
    psi = np.reshape(psi,[len(z),len(r)])
    plt.pcolor(r, z, psi)
    plt.plot(r_geom, z_geom, 'k-')
    plt.plot(r_left_target, z_left_target, 'g-')
    plt.plot(r_right_target, z_right_target, 'r-')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Magnetic Flux')
    # set the limits of the plot to the limits of the data
    plt.axis([r.min(), r.max(), z.min(), z.max()])
    plt.colorbar(label='Flux [Wb/rad')
    print( 'Saving psi function as psi.png ')
    plt.savefig('psi.png')
    plt.close()

    plt.contour(r,z,psi,100)
    plt.plot(r_geom, z_geom, 'k-')
    plt.plot(r_left_target, z_left_target, 'g-')
    plt.plot(r_right_target, z_right_target, 'r-')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Magnetic Flux Contours')
    # set the limits of the plot to the limits of the data
    plt.axis([r.min(), r.max(), z.min(), z.max()])
    plt.colorbar(label='Flux [Wb/rad]')
    print ('Saving psi contour as psiContour.png ')
    plt.savefig('psiContour.pdf')

    print('Take gradients of magnetic flux to produce magnetic field')
    [gradz,gradr] = np.gradient(np.array(psi),z[1]-z[0],r[1]-r[0])

    br = -gradz/r
    bz =gradr/r

    plt.close()
    plt.pcolor(r,z,br)
    plt.plot(r_geom, z_geom, 'k-')
    plt.plot(r_left_target, z_left_target, 'g-')
    plt.plot(r_right_target, z_right_target, 'r-')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Br')
    plt.colorbar(label='B field [T]')
    print ('Saving br profile as br.png ')
    plt.savefig('br.png')
    plt.close()

    plt.pcolor(r,z,bz)
    plt.plot(r_geom, z_geom, 'k-')
    plt.plot(r_left_target, z_left_target, 'g-')
    plt.plot(r_right_target, z_right_target, 'r-')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Bz')
    plt.colorbar(label='B field [T]')
    print ('Saving br profile as br.png ')
    plt.savefig('bz.png')

    Bp = np.sqrt(np.multiply(br,br) + np.multiply(bz,bz))
    bt = ml.repmat(btf*rtf/r,len(z),1)
    bt = np.reshape(bt,[len(z),len(r)])

    plt.close()
    plt.pcolor(r,z,bt)
    plt.plot(r_geom, z_geom, 'k-')
    plt.plot(r_left_target, z_left_target, 'g-')
    plt.plot(r_right_target, z_right_target, 'r-')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Bt')
    plt.colorbar(label='B field [T]')
    print( 'Saving bt profile as bt.png ')
    plt.savefig('bt.png')
    plt.close()

    rootgrp = netCDF4.Dataset("bField.nc", "w", format="NETCDF4")
    nrr = rootgrp.createDimension("nR", len(r))
    nzz = rootgrp.createDimension("nZ", len(z))
    brr = rootgrp.createVariable("br","f8",("nZ","nR"))
    btt = rootgrp.createVariable("bt","f8",("nZ","nR"))
    bzz = rootgrp.createVariable("bz","f8",("nZ","nR"))
    psii = rootgrp.createVariable("psi","f8",("nZ","nR"))
    rr = rootgrp.createVariable("r","f8",("nR"))
    zz = rootgrp.createVariable("z","f8",("nZ"))
    brr[:] = br
    btt[:] = bt
    bzz[:] = bz
    psii[:] = psi
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
def getBfield(rTarg,zTarg, \
              filename='/global/homes/t/tyounkin/atomIPS/atom-install-edison/solps-iter-data/Baseline2008-li0.70.x4.equ', \
              geometryFile='/global/homes/t/tyounkin/atomIPS/atom-install-edison/GITR/iter/iter_milestone/2d/input/iter2dRefinedOuterTarget.cfg',
              rmin=5.55,rmax=6.226,zmin=-4.6,zmax=-3.238):
    r,z,br, bz, bt,psi = readEquilibrium(filename)
    rSep,bAngle, bMag = interpolateBfield(r,z,br, bz, bt,psi,rTarg,zTarg,geometryFile,rmin,rmax,zmin,zmax)
    print( 'plotting bangle at target')
    plt.close()
    plt.plot(rSep,bAngle)
    plt.savefig('bangle.png')
    return rSep,bAngle,bMag

def process_solps_output_for_gitr(dakota_filename = '/Users/Alyssa/Dev/solps-iter-data/build/dakota', \
                                  nR = 500, nZ = 1000, plot_variables=0, \
                                  b2fstate_filename = '/Users/Alyssa/Dev/WEST/helium/b2fstate'):
    nIonSpecies, am, zamin, zn = get_solps_species(b2fstate_filename)

    dak = np.loadtxt(dakota_filename)
    dak = np.reshape(dak, (nR * nZ, -1))
    print('dak shape', dak.shape)

    rdak = np.unique(dak[:, 0])
    rdak = np.linspace(rdak[0],rdak[-1],len(rdak))
    zdak = np.unique(dak[:, 1])
    zdak = np.linspace(zdak[0],zdak[-1],len(zdak))

    te = get_dakota_variable(2,dak,rdak,zdak,nR,nZ,'te',plot_variables)
    off_grid_inds = np.where(te < -0.5)
    te[off_grid_inds] = 0.0;

    ne = get_dakota_variable(3,dak,rdak,zdak,nR,nZ,'ne',plot_variables)
    ne[off_grid_inds] = 0.0;

    ti = get_dakota_variable(4,dak,rdak,zdak,nR,nZ,'ti',plot_variables)
    ti[off_grid_inds] = 0.0;

    ni = np.zeros((nIonSpecies, nZ, nR))
    v_parallel = np.zeros((nIonSpecies, nZ, nR))

    for i in range(nIonSpecies):
        ni[i,:,:] = get_dakota_variable(5+i,dak,rdak,zdak,nR,nZ,'ni'+str(i),plot_variables)
        v_parallel[i,:,:] = get_dakota_variable(5+4*nIonSpecies+i,dak,rdak,zdak,nR,nZ,'v_parallel'+str(i),plot_variables)

    ni_total = np.zeros(( nZ, nR))
    v_parallel_total = np.zeros((nZ, nR))
    aveMass = np.zeros((nZ, nR))
    aveCharge = np.zeros((nZ, nR))

    for i in range(nIonSpecies):
        if zamin[i] > 0.0:
            ni_total = ni_total + np.reshape(ni[i, :, :], (nZ, nR))
            aveMass = aveMass + np.reshape(am[i] * ni[i, :, :], (nZ, nR))
            aveCharge = aveCharge + np.reshape(zamin[i] * ni[i, :, :], (nZ, nR))
            v_parallel_total = v_parallel_total + np.reshape(np.multiply(v_parallel[i, :, :], ni[i, :, :]), (nZ, nR))

    aveMass = np.divide(aveMass, ni_total)
    aveMass[off_grid_inds] = 0
    aveCharge = np.divide(aveCharge, ni_total)
    aveCharge[off_grid_inds] = 0
    ni_total[off_grid_inds] = 0
    v_parallel_total = np.divide(v_parallel_total, ni_total)
    v_parallel_total[off_grid_inds] = 0

    if plot_variables:
        plt.close()
        plt.pcolor(rdak, zdak, ni_total)
        plt.colorbar()
        plt.savefig('niTotal.png')
        plt.close()
        plt.pcolor(rdak, zdak, aveMass)
        plt.colorbar()
        plt.savefig('aveMass.png')
        plt.close()
        plt.pcolor(rdak, zdak, aveCharge)
        plt.colorbar()
        plt.savefig('aveCharge.png')

    br = get_dakota_variable(5+5*nIonSpecies,dak,rdak,zdak,nR,nZ,'br',plot_variables)
    bphi = get_dakota_variable(5+5*nIonSpecies+1,dak,rdak,zdak,nR,nZ,'bphi',plot_variables)
    bz = get_dakota_variable(5+5*nIonSpecies+2,dak,rdak,zdak,nR,nZ,'bz',plot_variables)

    print('br size',br.shape)
    vr,vt,vz = project_parallel_variable_xyz(v_parallel_total, br, bphi, bz,rdak,zdak, nR, nZ, 'v',plot_variables)
    vr[off_grid_inds] = 0.0;
    vt[off_grid_inds] = 0.0;
    vz[off_grid_inds] = 0.0;

    grad_ti = get_dakota_variable(5+ 5*nIonSpecies+4, dak, rdak, zdak, nR, nZ, 'grad_ti',plot_variables)
    grad_te = get_dakota_variable(5+ 5*nIonSpecies+5, dak, rdak, zdak, nR, nZ, 'grad_te',plot_variables)

    print('br size',br.shape)

    grad_ti_r,grad_ti_t,grad_ti_z = project_parallel_variable_xyz(grad_ti, br, bphi, bz,rdak,zdak, nR, nZ, 'grad_ti',plot_variables)
    grad_te_r,grad_te_t,grad_te_z = project_parallel_variable_xyz(grad_te, br, bphi, bz,rdak,zdak, nR, nZ, 'grad_te',plot_variables)
    grad_ti_r[off_grid_inds] = 0.0;
    grad_ti_t[off_grid_inds] = 0.0;
    grad_ti_z[off_grid_inds] = 0.0;
    grad_te_r[off_grid_inds] = 0.0;
    grad_te_t[off_grid_inds] = 0.0;
    grad_te_z[off_grid_inds] = 0.0;
    e_para = get_dakota_variable(5+ 5*nIonSpecies+6, dak, rdak, zdak, nR, nZ, 'e_para',plot_variables)
    e_perp = get_dakota_variable(5+ 5*nIonSpecies+7, dak, rdak, zdak, nR, nZ, 'e_perp',plot_variables)
    e_para[off_grid_inds] = 0.0;
    e_perp[off_grid_inds] = 0.0;
    profiles_filename = "profiles.nc"
    if os.path.exists(profiles_filename):
        os.remove(profiles_filename)

    rootgrp = netCDF4.Dataset(profiles_filename, "w", format="NETCDF4")
    nrr = rootgrp.createDimension("nR", len(rdak))
    nzz = rootgrp.createDimension("nZ", len(zdak))
    brr = rootgrp.createVariable("br", "f8", ("nZ", "nR"))
    btt = rootgrp.createVariable("bt", "f8", ("nZ", "nR"))
    bzz = rootgrp.createVariable("bz", "f8", ("nZ", "nR"))
    rr = rootgrp.createVariable("r", "f8", ("nR"))
    zz = rootgrp.createVariable("z", "f8", ("nZ"))
    brr[:] = br
    btt[:] = bphi
    bzz[:] = bz
    rr[:] = rdak
    zz[:] = zdak
    tee = rootgrp.createVariable("te", "f8", ("nZ", "nR"))
    nee = rootgrp.createVariable("ne", "f8", ("nZ", "nR"))
    tii = rootgrp.createVariable("ti", "f8", ("nZ", "nR"))
    nii = rootgrp.createVariable("ni", "f8", ("nZ", "nR"))
    nii1 = rootgrp.createVariable("ni1","f8", ("nZ", "nR"))
    nii2 = rootgrp.createVariable("ni2","f8", ("nZ", "nR"))
    mass = rootgrp.createVariable("mass", "f8", ("nZ", "nR"))
    charge = rootgrp.createVariable("charge", "f8", ("nZ", "nR"))
    v_para1 = rootgrp.createVariable("v_para1", "f8", ("nZ", "nR"))
    v_para2 = rootgrp.createVariable("v_para2", "f8", ("nZ", "nR"))
    vrr = rootgrp.createVariable("vr", "f8", ("nZ", "nR"))
    vzz = rootgrp.createVariable("vz", "f8", ("nZ", "nR"))
    vpp = rootgrp.createVariable("vp", "f8", ("nZ", "nR"))
    Err = rootgrp.createVariable("Er", "f8", ("nZ", "nR"))
    Ett = rootgrp.createVariable("Et", "f8", ("nZ", "nR"))
    Ezz = rootgrp.createVariable("Ez", "f8", ("nZ", "nR"))
    teer = rootgrp.createVariable("gradTeR", "f8", ("nZ", "nR"))
    teez = rootgrp.createVariable("gradTeZ", "f8", ("nZ", "nR"))
    teey = rootgrp.createVariable("gradTeY", "f8", ("nZ", "nR"))
    tiir = rootgrp.createVariable("gradTiR", "f8", ("nZ", "nR"))
    tiiz = rootgrp.createVariable("gradTiZ", "f8", ("nZ", "nR"))
    tiiy = rootgrp.createVariable("gradTiY", "f8", ("nZ", "nR"))
    tee[:] = te
    nee[:] = ne
    tii[:] = ti
    nii[:] = ni_total
    nii1[:] = ni[1]
    nii2[:] = ni[2]
    mass[:] = aveMass
    charge[:] = aveCharge
    v_para1[:] = v_parallel[1]
    v_para2[:] = v_parallel[2]
    vrr[:] = vr
    vpp[:] = vt
    vzz[:] = vz
    Err[:] = e_para
    Ett[:] = 0 * e_para
    Ezz[:] = e_perp
    teer[:] = grad_te_r
    teey[:] = grad_te_t
    teez[:] = grad_te_z
    tiir[:] = grad_ti_r
    tiiy[:] = grad_ti_t
    tiiz[:] = grad_ti_z
    rootgrp.close()

def project_parallel_variable_xyz(v_parallel_total,br,bphi,bz,rdak,zdak,nR,nZ,title='v',plot_variables=0):
    b_total = np.sqrt(np.multiply(br,br) + np.multiply(bphi,bphi) + np.multiply(bz,bz))

    vr = np.divide(np.multiply(br,v_parallel_total),b_total)
    vt = np.divide(np.multiply(bphi,v_parallel_total),b_total)
    vz = np.divide(np.multiply(bz,v_parallel_total),b_total)

    off_grid_inds = np.where(br == -1)
    if plot_variables:
        vr[off_grid_inds] = float("nan")
        vt[off_grid_inds] = float("nan")
        vz[off_grid_inds] = float("nan")

        plt.close()
        plt.pcolor(rdak, zdak, vr)
        plt.colorbar()
        plt.savefig(title+'_r.png')
        plt.close()
        plt.pcolor(rdak, zdak, vt)
        plt.colorbar()
        plt.savefig(title+'_t.png')
        plt.close()
        plt.pcolor(rdak, zdak, vz)
        plt.colorbar()
        plt.savefig(title+'_z.png')

        vr[off_grid_inds] = -1
        vt[off_grid_inds] = -1
        vz[off_grid_inds] = -1

    return vr,vt,vz

def get_dakota_variable(index,dak,rdak,zdak,nR,nZ,title='title',plot_variables=0):
    variable = np.reshape(dak[:, index], (nZ, nR))

    if plot_variables:
        off_grid_inds = np.where(variable == -1)
        variable[off_grid_inds] = float("nan");
        plt.close()
        plt.pcolor(rdak, zdak, np.reshape(variable, (nZ, nR)))
        plt.colorbar(orientation='vertical')
        plt.savefig(title+'.png')
        plt.close()
        variable[off_grid_inds] = -1;

    return variable

def intersection(a_x, a_y, b_x, b_y):
    i_a = np.array(0)
    i_b = np.array(0)

    for i in range(len(b_x)):
        for j in range(len(a_x)):
            if (b_x[i] == a_x[j] and b_y[i] == a_y[j]):
                i_a = np.append(i_a, j)
                i_b = np.append(i_b, i)

    i_a = np.delete(i_a, 0)
    i_b = np.delete(i_b, 0)

    return i_a, i_b

def find_strike_points(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry'):
    nx, ny, crx, cry, region = read_b2f_geometry(solps_geometry_filename)

    print('crx_size',crx.shape)
    print('region_size',region.shape)
    bottom_left = 0
    bottom_right = 1

    region_3_position = bottom_right
    region_4_position = bottom_left

    crx_region_3 = crx[:,:,region_3_position]
    cry_region_3 = cry[:,:,region_3_position]

    crx_region_4 = crx[:,:,region_4_position]
    cry_region_4 = cry[:,:,region_4_position]
    print('crx_region_3 shape',crx_region_3.shape)
    print('crx_region_4 shape',crx_region_4.shape)

    condition_region_3 = region == 3
    condition_region_4 = region == 4

    region_3_indx = np.where(condition_region_3)
    region_4_indx = np.where(condition_region_4)

    print('region_3_indx shape',region_3_indx)
    x_region_3 = crx_region_3[region_3_indx]
    y_region_3 = cry_region_3[region_3_indx]
    x_region_4 = crx_region_4[region_4_indx]
    y_region_4 = cry_region_4[region_4_indx]

    plt.close()
    plt.scatter(x_region_3,y_region_3)
    plt.scatter(x_region_3,y_region_3)
    plt.savefig('region34.png')
    plt.close()

    i_a,i_b = intersection(x_region_3,y_region_3,x_region_4,y_region_4)

    index = np.argmax(y_region_3[i_a])
    print('index',index)
    print('ia_index',i_a[index])
    row_col = np.unravel_index(i_a[index],crx_region_3.shape)

    x_x_point = crx_region_3[row_col[0],row_col[1]]
    y_x_point = cry_region_3[row_col[0],row_col[1]]

    x_inner_strikepoint = crx_region_3[0,row_col[1]]
    y_inner_strikepoint = cry_region_3[0,row_col[1]]

    x_outer_strikepoint = crx_region_3[-1,row_col[1]]
    y_outer_strikepoint = cry_region_3[-1,row_col[1]]

    topcut = read_b2f_variable(solps_geometry_filename,'topcut')
    topcut = int(topcut)+1; # to go from solps index (-1) to python index (0)

    print('topcut',topcut)
    x_inner_strikepoint = crx[1,topcut,bottom_left]
    y_inner_strikepoint = cry[1,topcut,bottom_left]
    x_outer_strikepoint = crx[-1,topcut,bottom_left]
    y_outer_strikepoint = cry[-1,topcut,bottom_left]

    print('sp xy',x_outer_strikepoint,y_outer_strikepoint)
    print('sp xy crx',crx[-1,topcut,:])
    print('sp xy cry',cry[-1,topcut,:])
    return x_x_point,y_x_point, \
           x_inner_strikepoint ,y_inner_strikepoint, \
           x_outer_strikepoint ,y_outer_strikepoint, topcut



def get_target_coordinates(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry'):

    # Get number of mesh elements in x and y (SOLPS coordinates), nx, ny.
    # As well as the coordinates of the corners, crx, cry, 
    # and the region number from solps_geometry_filename
    nx, ny, crx, cry, region = read_b2f_geometry(solps_geometry_filename)
    print('Number of SOLPS gridpoints in x: ', nx)
    print('Number of SOLPS gridpoints in y: ', ny)

    geom_shape = crx.shape
    bottom_left = 0
    bottom_right = 1
    top_left = 2;
    top_right = 3;

    r_inner_target = crx[0,1:,bottom_right]
    z_inner_target = cry[0,1:,bottom_right]
    r_outer_target = crx[-1,1:,bottom_left]
    z_outer_target = cry[-1,1:,bottom_left]
    print('Number of inner and outer target points (should be ny-1): ',r_inner_target.size)
    # Therefore there should be ny-2 line segments from the solps mesh to
    # be introduced to the gitr geometry

    return r_inner_target,z_inner_target, \
           r_outer_target,z_outer_target


def read_b2f_geometry(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry'):
    nxny = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='nx,ny')
    nx = int(nxny[0]+2)
    ny = int(nxny[1]+2)
    #print('nx,ny',nx,ny)
    crx = np.zeros((nx, ny, 4))
    cry = np.zeros((nx, ny, 4))

    crx_long = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='crx')
    cry_long = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='cry')

    for i in range(4):
        crx[:,:,i] = np.transpose(np.reshape(crx_long[i*nx*ny:(i+1)*nx*ny],(ny,nx)))
        cry[:,:,i] = np.transpose(np.reshape(cry_long[i*nx*ny:(i+1)*nx*ny],(ny,nx)))

    #print('crx shape',crx.shape)
    region = read_b2f_variable(solps_geometry_filename, \
                            field_name='region')
    #print('firstwhatever',region[0:nx*ny])
    region = np.transpose(np.reshape(region[0:nx*ny],(ny,nx)))

    return nx,ny,crx,cry,region
def read_b2f_variable(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry', \
                      field_name = 'crx'):
    f = open(solps_geometry_filename, 'r')
    txt = f.readlines()
    f.close()

    field_start = 0
    field_end = 0
    found = 0

    for count, line in enumerate(txt):
        if found == 0:
            if '*cf' in line:
                words = line.split()
                if words[-1] == field_name:
                    field_start = count+1
                    found = 1;
        elif found == 1:
            if '*cf' in line:
                field_end = count
                found = 2
        elif found == 2:
            break

    field = [];
    txt_list = txt[field_start:field_end]
    for sublist in txt_list:
        split_sublist = sublist.split()
        for element in split_sublist:
            field.append(element)

    field = np.array(field)
    field = field.astype(np.float)

    return field

def get_solps_species(solps_state_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fstate'):
    f = open(solps_state_filename, 'r')
    txt = f.readlines()[:200]
    f.close()

    znLine = 0
    zn = ''
    for count, line in enumerate(txt):
        if znLine:
            if '*' in line:
                break
            else:
                zn = ''.join((zn, line))
        if 'zn' in line:
            znLine = count

    zaminLine = 0
    zamin = ''
    for count, line in enumerate(txt):
        if zaminLine:
            if '*' in line:
                break
            else:
                zamin = ''.join((zamin, line))
        if 'zamin' in line:
            zaminLine = count

    amLine = 0
    am = ''
    for count, line in enumerate(txt):
        if amLine:
            if '*' in line:
                break
            else:
                am = ''.join((am, line))
        if 'am ' in line:
            amLine = count

    zn = zn.split()
    zn = [float(i) for i in zn]

    zamin = zamin.split()
    zamin = [float(i) for i in zamin]

    am = am.split()
    am = [float(i) for i in am]
    # Section to manually add in Tritium
    #zn.insert(2, 1.0)
    #zamin.insert(2, 0.0)
    #am.insert(2, 3.0)
    #zn.insert(3, 1.0)
    #zamin.insert(3, 1.0)
    #am.insert(3, 3.0)
    nIonSpecies = len(zn)
    species_index = [int(i) for i in range(nIonSpecies)]

    array = np.array((species_index,zn,am,zamin))
    array = np.transpose(array)
    print(array)
    np.savetxt('speciesList.txt',array,header='SpeciesIndex   Z   Mass   Charge', delimiter=' ') 
    #print('nIonSpecies', nIonSpecies)
    #species_file = open("speciesList.txt", "w")
    #species_file.write('SpeciesIndex   Z   Mass   Charge\n')
    #print('Existing species for current SOLPS run:\n')
    #print('SpeciesIndex   Z   Mass   Charge\n')
    #for i in range(nIonSpecies):
    #    print('%f %f %f %f \n' % (i, zn[i], am[i], zamin[i]))
    #    species_file.write('%f %f %f %f \n' % (i, zn[i], am[i], zamin[i]))

    #species_file.close()

    return nIonSpecies, am,zamin, zn

def read_target_file(filename = '/Users/Alyssa/Dev/solps-iter-data/build/rightTargOutput' ):
    # target files contain
    # r,z,ti,
    # ni for nSpecies,
    # flux ro nSpeces,
    # te, ne
    target = np.loadtxt(filename)
    target_shape = target.shape
    nSpecies = int((target_shape[1] - 5)/2)
    print('nSpecies',nSpecies)
    r = target[:, 0]
    z = target[:, 1]
    ti = target[:,2]
    ni = np.zeros((target_shape[0],nSpecies))
    flux = np.zeros((target_shape[0],nSpecies))
    for i in range(nSpecies):
        ni[:,i] = target[:,3+i]
        flux[:,i] = target[:,3+nSpecies+i]

    te = target[:,3+2*nSpecies]
    ne = target[:,3+2*nSpecies+1]

    return r,z,ti,ni,flux,te,ne

def make_solps_targ_coord_file(gitr_geom_filename='/Users/Alyssa/Dev/GITR/west/helium/input/gitrGeometry.cfg', \
    solps_geom = '/Users/Alyssa/Dev/WEST/baserun/b2fgmtry', \
    coords_file = '/Users/Alyssa/Dev/GITR/west/helium/input/right_target_coordinates.txt', \
    right_target_filename= '/Users/Alyssa/Dev/solps-iter-data/build/rightTargOutput'):
    r, z, ti, ni, flux, te, ne = read_target_file(right_target_filename)

    x_x_point, y_x_point, \
    x_inner_strikepoint, y_inner_strikepoint, \
    x_outer_strikepoint, y_outer_strikepoint, topcut = find_strike_points(solps_geom)
    topcut = topcut - 1; #Because the target coordinates already strip off the ghost cell

    r_left_target,z_left_target,r_right_target,z_right_target = get_target_coordinates(solps_geom)
    print('r strikepoint and target coords',x_outer_strikepoint, r_right_target)
    print('r strikepoint length',len(r_right_target))

    r_minus_r_sep = 0*r_right_target

    lengths = np.sqrt(np.multiply((r_right_target[1:] -r_right_target[0:-1]),(r_right_target[1:] -r_right_target[0:-1])) + \
                      np.multiply((z_right_target[1:] -z_right_target[0:-1]),(z_right_target[1:] -z_right_target[0:-1])) )

    for i in range(topcut+1,len(r_right_target)):
        r_minus_r_sep[i] = r_minus_r_sep[i-1]+lengths[i-1]

    for i in range(topcut-1,-1,-1):
        r_minus_r_sep[i] = r_minus_r_sep[i+1]-lengths[i]

    print('r_minus_r_sep',r_minus_r_sep)
    #tol = 1e-4
    #condition = np.abs(r - x_outer_strikepoint) <= tol
    #condition2 = np.abs(z - y_outer_strikepoint) <= tol

    #rz_ind = np.where(condition & condition2)
    #print('strikepoint',x_outer_strikepoint,y_outer_strikepoint)
    #print('strikepoint rzind',rz_ind)
    #print('strikepoint rz',r,z)
    plt.close()
    plt.plot(r,z)
    plt.scatter(r,z)
    plt.scatter(r_right_target,z_right_target)
    #plt.scatter(x_outer_strikepoint,y_outer_strikepoint)
    #plt.axis([5.5,5.6,-4.42, -4.36])
    plt.savefig('rz_targ.png')
    plt.show()

    with io.open(gitr_geom_filename) as f:
        config = libconf.load(f)

    x1 = np.array(config.geom.x1)
    x2 = np.array(config.geom.x2)
    z1 = np.array(config.geom.z1)
    z2 = np.array(config.geom.z2)

    i_a, i_b = intersection(x1, z1, r_right_target, z_right_target)
    print('i_a',i_a)
    A = np.zeros((len(i_a),10))
    A[:,0] = r_right_target
    A[:,1] = z_right_target
    A[:,2] = r_minus_r_sep
    A[:,3] = i_a
    A[:,4] = x1[i_a]
    A[:,5] = x2[i_a]
    A[:,6] = z1[i_a]
    A[:,7] = z2[i_a]
    A[:,8] = 0.5*(x1[i_a]+x2[i_a])
    A[:,9] = 0.5*(z1[i_a]+z2[i_a])

    np.savetxt(coords_file,A,header='r,z,r_minus_r_sep,gitr_index,x1,x2,z1,z2,xmid,zmid')

def make_solps_targ_file(solps_geom = '/Users/Alyssa/Dev/WEST/baserun/b2fgmtry', \
    b_field_file = '/Users/Alyssa/Dev/WEST/baserun/west_54034_10p2s_mag.X4.equ', \
    coords_file = '/Users/Alyssa/Dev/GITR/west/helium/input/right_target_coordinates.txt', \
    right_target_filename= '/Users/Alyssa/Dev/solps-iter-data/build/rightTargOutput'):

    r, z, ti, ni, flux, te, ne = read_target_file(right_target_filename)

    r_left_target,z_left_target,r_right_target,z_right_target = get_target_coordinates(solps_geom)

    rr, zz, r_minus_r_sep, gitr_ind, r_midpoint, z_midpoint = read_targ_coordinates_file(coords_file)

    slope = np.zeros(len(rr)-1)
    rise = zz[1:] - zz[0:-1]
    run = rr[1:] - rr[0:-1]
    slope = np.divide(rise,run)
    perp_slope = -np.sign(slope)/np.absolute(slope)
    rPerp = np.divide(1.0,np.sqrt(np.multiply(perp_slope,perp_slope)+1))
    zPerp = 1.0*np.sqrt(1-np.multiply(rPerp,rPerp))
    length = np.sqrt(np.multiply(rise,rise)+np.multiply(run,run))
    #print('slope',slope)
    #print('lenslope',len(slope))
    #z_midpoint = 0.5*(z_no_guard_cells[1:] + z_no_guard_cells[0:-1])
    #r_midpoint = 0.5*(r_no_guard_cells[1:] + r_no_guard_cells[0:-1])
    print('rmid', r_midpoint)
    print('zmid', z_midpoint)

    r, z, br, bz, bt, psi = readEquilibrium(b_field_file)
    
    btot = np.sqrt(np.multiply(br,br) + np.multiply(bz,bz) + np.multiply(bt,bt))

    #f = scii.interp2d(r,z,btot)
    #btarg = f(r_midpoint,z_midpoint)
    grid_r, grid_z = np.meshgrid(r,z)
    print(grid_r.shape, grid_z.shape, btot.shape)
    btarg = scii.griddata((grid_r.flatten(),grid_z.flatten()), btot.flatten(), (r_midpoint, z_midpoint), method='linear')
    brtarg = scii.griddata((grid_r.flatten(),grid_z.flatten()), br.flatten(), (r_midpoint, z_midpoint), method='linear')
    bztarg = scii.griddata((grid_r.flatten(),grid_z.flatten()), bz.flatten(), (r_midpoint, z_midpoint), method='linear')
    
    angle = 180.0/np.pi*np.arccos(np.divide(np.multiply(brtarg,rPerp) + np.multiply(bztarg,zPerp),btarg))

    len_rmid = len(r_midpoint)-2
    A = np.zeros((len(r_midpoint)-2,3))
    A[:,0] = r_minus_r_sep[1:-2] + 0.5*length[1:-1]
    A[:,1] = r_midpoint[1:-1]
    A[:,2] = z_midpoint[1:-1]

    A = np.hstack((A,np.reshape(te,(len_rmid,1))))
    A = np.hstack((A,np.reshape(ti,(len_rmid,1))))
    A = np.hstack((A,flux))
    A = np.hstack((A,ni))
    A = np.hstack((A,np.reshape(btarg[1:-1],(len_rmid,1))))
    A = np.hstack((A,np.reshape(angle[1:-1],(len_rmid,1))))

    print('fshpe',flux.shape)
    np.savetxt('solpsTarg.txt',A,delimiter=',',header='R-Rsep, r, z, Te, Ti, Flux (for each species), n (for each species), Btot, Bangle')

def read_targ_coordinates_file(filename = 'right_target_coordinates.txt'):
    
    data = np.loadtxt(filename, skiprows=0)

    r_right_target = data[:,0]
    #r_right_target = np.append(r_right_target,data[-1,4])
    #print('rrt length',len(r_right_target),data[-1,4])
    z_right_target = data[:,1]
    #z_right_target = np.append(z_right_target,data[-1,6])
    r_minus_r_sep = data[:,2]
    #len_last_segment = np.sqrt((r_right_target[-1] - r_right_target[-2])**2 + (z_right_target[-1] - z_right_target[-2])**2)
    #print('rmrs length',len(r_minus_r_sep),len_last_segment)
    #print('data length',len(data[:,0]),len(data[:,1]),len(data[:,2]))
    #r_minus_r_sep = np.append(data[:,2], len_last_segment)
    
    gitr_ind = data[:,3]
    r_midpoint = data[0:-1,8]
    z_midpoint = data[0:-1,9]

    return r_right_target, z_right_target, r_minus_r_sep, gitr_ind, r_midpoint, z_midpoint



if __name__ == "__main__":   
    #rTarg = np.linspace(5,6.5,100)
    #zTarg=np.linspace(0,1,100)
    #getBfield(rTarg,zTarg,"/Users/tyounkin/Dissertation/ITER/mq3/final/Baseline2008-li0.70.x4.equ","/Users/tyounkin/Code/gitr2/iter/iter_milestone/2d/input/iterGeom2DdirBe0.cfg")
    process_solps_output_for_gitr()
    #get_solps_species()
    readEquilibrium()
    #read_b2f_geometry()
    #find_strike_points()
    #read_target_file()
    #get_target_coordinates()
    #make_solps_targ_coord_file()
    #make_solps_targ_file()
    #make_solps_targ_file(gitr_geom_filename='gitr_geometry.cfg', \
    #solps_geom = '/project/projectdirs/m1709/psi-install-cori/solps_data/mq3/b2fgmtry', \
    #right_target_filename= 'rightTargOutput')
    #make_solps_targ_file_txt(solps_geom='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry',b_field_file = '/Users/tyounkin/Dissertation/ITER/mq3/solps/Baseline2008-li0.70.x4.equ')
    #solps_geom = '/Users/tyounkin/postDoc/DOE-West/Deuterium/WEST_D_run1/baserun/b2fgmtry', \
    #b_field_file = '/Users/tyounkin/postDoc/DOE-West/Deuterium/WEST_D_run1/baserun/west_54034_10p2s_mag.X4.equ', \
    #coords_file = 'right_target_coordinates.txt', \
    #right_target_filename= 'rightTargOutput')
