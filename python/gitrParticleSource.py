from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
import gitr
import scipy.interpolate as scii
import netCDF4

def particleSource(spylFile = 'SpylPerSpecies.dat',solps_path='solpsTarg.txt'):
    spyl = np.loadtxt(spylFile)
    print 'spyl',spyl.shape
    print 'spyl',spyl
    SOLPS = np.loadtxt(solps_path, dtype='float',skiprows=1,delimiter=' ')
    print 'solps',SOLPS.shape
    spFlux = 0*spyl
    rmrs = SOLPS[:,0]
    r = SOLPS[:,1]
    z = SOLPS[:,2]
    #for i in range(len(rmrs)):
    fi  = abs(SOLPS[:,[21+5,21+5,28,29,31,32,33,34,36,37,38,39,40,41,42,43,44,45]]);
	#spFlux[i,:] = np.multiply(abs(fi),spyl[i,:])
    spFlux = np.reshape(np.sum(fi,axis=1),(36,1))*spyl
    plt.close()
    plt.plot(rmrs,spFlux[:,0:8])
    plt.plot(rmrs,spFlux[:,8:16], linestyle='dashed')
    plt.plot(rmrs,spFlux[:,16:19], linestyle='dotted')
    plt.title('Eroded Flux By Background Species')
    plt.ylabel('W Flux [m-2s-1]')
    plt.xlabel('R-Rsep [m]')
    plt.legend(['D_1+','T_1+','He_1+','He_2+','Be_1+','Be_2+','Be_3+','Be_4+','Ne_1+','Ne_2+','Ne_3+','Ne_4+','Ne_5+','Ne_6+','Ne_7+','Ne_8+','Ne_9+','Ne_10+'],loc=1)
    #plt.savefig('totalSputt.png')
    plt.show()
    heTotal = np.sum(spFlux[:,2:4],axis=1)
    beTotal = np.sum(spFlux[:,4:8],axis=1)
    neTotal = np.sum(spFlux[:,8:18],axis=1)
    plt.close()
    plt.plot(rmrs,spFlux[:,0])
    plt.plot(rmrs,spFlux[:,1])
    plt.plot(rmrs,heTotal)
    plt.plot(rmrs,beTotal)
    plt.plot(rmrs,neTotal)
    plt.title('Eroded Flux By Background Species')
    plt.ylabel('W Flux [m-2s-1]')
    plt.xlabel('R-Rsep [m]')
    plt.legend(['D','T','He','Be','Ne'])
    plt.savefig('totalSputtFlux2.png')
    totalFlux = spFlux[:,0]+spFlux[:,1]+heTotal+neTotal
    print('Total flux',totalFlux)
    #print('Total flux',totalFlux.shape())
    thisRmrs = scii.interpn([z],rmrs, [-4.3])
    localFlux = scii.interpn([rmrs],totalFlux,[thisRmrs])
    print('Local flux',localFlux)
    localFlux = scii.interpn([z],totalFlux,[-4.3])
    print('Local flux',localFlux)
    x1,x2,x3,y1,y2,y3,z1,z2,z3,area,Z,surfIndArray,surf,a,b,c,d,plane_norm=gitr.read3dGeom(filename='/Users/tyounkin/Code/gitr2/iter/iter_milestone/3d/input/iterRefinedTest.cfg')
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
    print "number of w surfaces",len(particleSurfs[0])
    print "number of particle surfaces",len(erodedSurfs)
    particleCDF = np.cumsum(surfaceParticleRate[particleSurfs[0]])
    particleCDF = 1/particleCDF[-1]*particleCDF
    print('particleCDF',particleCDF)
    print('particleCDF len',len(particleCDF))
    nParticles = 10000
    particleCDFmat = np.tile(particleCDF,(nParticles,1))
    print('particleCDFmat',particleCDFmat.shape)
    rand1 = np.random.rand(nParticles)
    diff = particleCDFmat.T - rand1
    diff[diff<0.0] = 100
    mins = np.argmin(diff,axis=0)
    print('mins',mins)
    xx = []
    yy=[]
    zz=[]
    vx = []
    vy=[]
    vz=[]
    buff=1e-6
    pE = 4
    pVel = np.sqrt(2*pE*1.602e-19/184/1.66e-27)
    print('particle velocity is ',pVel)

    for i in range(nParticles):
        xyz = sampleTriangle(x = np.array([x1[erodedSurfs[mins[i]]],x2[erodedSurfs[mins[i]]],x3[erodedSurfs[mins[i]]]]),
	y = np.array([y1[erodedSurfs[mins[i]]],y2[erodedSurfs[mins[i]]],y3[erodedSurfs[mins[i]]]]),
	z = np.array([z1[erodedSurfs[mins[i]]],z2[erodedSurfs[mins[i]]],z3[erodedSurfs[mins[i]]]]))
	xx.append(xyz[0] - a[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*buff)
	yy.append(xyz[1] - b[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*buff)
	zz.append(xyz[2] - c[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*buff)
	vx.append( - a[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*pVel)
	vy.append( - b[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*pVel)
	vz.append( - c[erodedSurfs[mins[i]]]/plane_norm[erodedSurfs[mins[i]]]*pVel)
        
    #xx = x1[erodedSurfs[mins]]
    #yy = y1[erodedSurfs[mins]]
    #zz = z1[erodedSurfs[mins]]
    print('xx',xx)
    print('yy',yy)
    print('zz',zz)

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
    print 'mininum z',z.min()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(xc, yc, zc, c=surfaceParticleRate, marker='o')
    ax.scatter(xx, yy, zz, c='r', marker='o')
    #ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
    plt.show()
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
if __name__ == "__main__":
    particleSource()
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
