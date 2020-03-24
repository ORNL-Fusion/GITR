# python library tools for gitr

from distutils.dir_util import copy_tree
import sys
# sys.path.append('/home/tqd/code/netcdf4-python')
import netCDF4
import numpy as np
# import Tkinter
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
# matplotlib.use('agg')
# import cv2
import io, libconf
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pylab as pl
import scipy as sp
import math
import os
import shutil
import math
import scipy.interpolate as scii


def copy_folder(from_folder, to_folder):
    copy_tree(from_folder, to_folder)


def nc_show(filename):
    ncFile = netCDF4.Dataset(filename, "r")
    print("File : ", filename, " format ", ncFile.file_format)
    # print("nLines ", surfFile.dimensions.keys())
    print("DIMENSIONS: ")
    for key in ncFile.dimensions:
        print(key, ncFile.dimensions[key].size)

    print("VARIABLES: ")
    for key in ncFile.variables:
        print(key, ncFile.variables[key].dimensions, ncFile.variables[key][:])


def depositedEdist():
    ncFile = netCDF4.Dataset("surface.nc", "r")
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
    totalEdist = np.sum(Edist[surfaceIndices][:][:], axis=0)
    eadistSurf1 = Edist[0][:][:]
    edistSurf1 = np.sum(eadistSurf1, axis=0)
    print('total Edist', totalEdist.shape)
    EdistOnly = np.sum(totalEdist, axis=0)
    AdistOnly = np.sum(totalEdist, axis=1)
    totalImpacts = np.sum(EdistOnly, axis=0)
    print('total Impacts ', totalImpacts)
    E = np.linspace(0, 1000, 200)
    plt.figure(1, figsize=(10, 6), dpi=80)
    # plt.plot(E,edistSurf1)
    plt.plot(gridE, EdistOnly)
    # plt.show()
    plt.savefig('image1.png')
    plt.figure(2, figsize=(10, 6), dpi=80)
    plt.plot(gridA, AdistOnly)
    plt.savefig('image2.png')
    # image = cv2.imread("image1.png")
    # cv2.imshow("Image", image)
    for i in range(1, len(gridA)):
        writeEDfile("angle" + str(i) + ".ED1", gridE, EdistOnly)

    writeEDfile("angularDist.dat", gridA, AdistOnly)


def writeEDfile(name, gridE, Edist):
    if (gridE[0] == 0):
        dE = gridE[1] - gridE[0]
        gridE = gridE + 0.5 * dE
    datafile_id = open(name, 'w+')
    data = np.array([gridE, Edist])
    data = data.T
    np.savetxt(datafile_id, data, fmt=['%.4f', '%.4f'])
    # here the ascii file is populated.
    datafile_id.close()


def modifyInputParam(filename="input/gitrInput.cfg", nT=100, nP=1000):
    with io.open(filename) as f:
        config = libconf.load(f)
    config['timeStep']['nT'] = nT
    config['impurityParticleSource']['nP'] = nP
    with io.open(filename, 'w') as f:
        libconf.dump(config, f)


def plot2dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)

    x1 = np.array(config.geom.x1)
    x2 = np.array(config.geom.x2)
    z1 = np.array(config.geom.z1)
    z2 = np.array(config.geom.z2)
    Z = np.array(config.geom.Z)
    length = np.array(config.geom.length)
    y1 = config.geom.y1
    y2 = config.geom.y2
    ys1 = np.ones(x1.size) * y1
    ys2 = np.ones(x1.size) * y2
    # if plt.fignum_exists(num=1):
    plt.plot(np.append(x1, x1[0]), np.append(z1, z1[0]), linewidth=2.0, color='k')
    # else:
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
    return x1, x2, z1, z2, length, Z


def read3dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)

    x1 = np.array(config.geom.x1)
    x2 = np.array(config.geom.x2)
    x3 = np.array(config.geom.x3)
    y1 = np.array(config.geom.y1)
    y2 = np.array(config.geom.y2)
    y3 = np.array(config.geom.y3)
    z1 = np.array(config.geom.z1)
    z2 = np.array(config.geom.z2)
    z3 = np.array(config.geom.z3)
    area = np.array(config.geom.area)
    surf = np.array(config.geom.surface)
    Z = np.array(config.geom.Z)
    materialSurfaceInidces = np.nonzero(Z)
    surfIndArray = np.asarray(materialSurfaceInidces)
    a = np.array(config.geom.a)
    b = np.array(config.geom.b)
    c = np.array(config.geom.c)
    d = np.array(config.geom.d)
    plane_norm = np.array(config.geom.plane_norm)
    print('Number of W surfaces ', surfIndArray.size)
    return x1, x2, x3, y1, y2, y3, z1, z2, z3, area, Z, surfIndArray, surf, a, b, c, d, plane_norm


def iter2dProcessing():
    x1, x2, z1, z2, length, Z = plot2dGeom('input/iterRefinedTest.cfg')
    plt.close()
    grossDep, grossEro, sumWeightStrike, E, A, EAdist, surfaceNumbers, surfEdist = nc_readSurface()
    plt.close()
    netErosion = grossDep - grossEro
    colormap = plt.cm.bwr
    normalize = matplotlib.colors.Normalize(vmin=-1., vmax=1.)
    plt.subplot(2, 1, 1)
    for i in range(0, len(surfaceNumbers)):
        plt.plot([x1[surfaceNumbers[i]], x2[surfaceNumbers[i]]], [z1[surfaceNumbers[i]], z2[surfaceNumbers[i]]],
                 color=colormap(0.))
        plt.hold(True)
    plt.subplot(2, 1, 2)
    plt.scatter(surfaceNumbers, netErosion)
    plt.savefig('surfaces.png')
    plt.close()
    x, y, r, z, charge = nc_plotPositions('output/positions.nc')
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
        thisSurf = surfaceNumbers[i + 4]
        lastSurf = surfaceNumbers[i + 3]
        print(thisSurf)
        print(z2[thisSurf])
        if (i == 0):
            rMinusRsep[i] = z2[thisSurf] - zSep
        else:
            rMinusRsep[i] = rMinusRsep[i - 1] + length[lastSurf]
    plt.close('all')
    s1 = plt.scatter(rMinusRsep, grossEro[range(4, 164)], c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep, grossDep[range(4, 164)], c='red')
    s3 = plt.scatter(rMinusRsep, netErosion[range(4, 164)], c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosion.png')
    plt.close('all')
    s1 = plt.scatter(rMinusRsep, np.log10(grossEro[range(4, 164)]), c='blue')
    plt.hold(True)
    s2 = plt.scatter(rMinusRsep, np.log10(grossDep[range(4, 164)]), c='red')
    s3 = plt.scatter(rMinusRsep, np.sign(netErosion[range(4, 164)]) * np.log10(np.absolute(netErosion[range(4, 164)])),
                     c='green')
    plt.hold(True)
    plt.legend((s1, s2, s3), ('Gross Erosion', 'Gross Deposition', 'netErosion'))
    plt.savefig('targetErosionlog.png')
    plt.close()


def printHeDist(path='', z=[-4.1, -4.0]):
    ncFile = netCDF4.Dataset("input/iterHeDist.nc", "r")
    nPoints = len(ncFile.dimensions['nPoints'])
    nE = len(ncFile.dimensions['nE'])
    nA = len(ncFile.dimensions['nA'])

    zPoints = np.array(ncFile.variables['z'])
    flux = np.array(ncFile.variables['flux'])
    EAdist = np.array(ncFile.variables['EAdist'])
    gridE = np.array(ncFile.variables['gridE'])
    gridA = np.array(ncFile.variables['gridA'])
    heFlux = np.zeros(len(z))
    os.mkdir('GITRoutput')
    for i in range(len(z)):
        idx = (np.abs(zPoints - z[i])).argmin()
        print('index', idx)
        heFlux[i] = flux[idx]
        thisDist = EAdist[:, :, idx]  # first dimension is angle, second is energy
        thisDist = np.transpose(thisDist)
        print('thisDistHe size', thisDist.shape)
        # print(thisDist)
        Aweight = np.sum(thisDist, axis=0)
        gitrDir = 'GITRoutput/' + 'gitr' + str(i)
        os.mkdir(gitrDir)
        np.savetxt(gitrDir + '/gitrFluxE.dat', gridE)
        np.savetxt(gitrDir + '/gitrFluxAweight.dat', Aweight)
        np.savetxt(gitrDir + '/gitrFluxA.dat', gridA[:-1])
        np.savetxt(gitrDir + '/gitrFluxEAdist.dat', thisDist)

        if (path != ''):
            np.savetxt(path + '/' + gitrDir + '/gitrFluxE.dat', gridE)
            np.savetxt(path + '/' + gitrDir + '/gitrFluxAweight.dat', Aweight)
            np.savetxt(path + '/' + gitrDir + '/gitrFluxA.dat', gridA)
            np.savetxt(path + '/' + gitrDir + '/gitrFluxEAdist.dat', thisDist)

        for j in range(0, nA):
            # print('Evec size ', gridE.size)
            # print('EAdist[:,i] size ', EAdist[:,i].size)
            edOut = np.column_stack((gridE, thisDist[:, j]))
            np.savetxt(gitrDir + '/dist' + str(j) + '.dat', edOut)
    return gitrDir, heFlux


def iter3dProcessing(path='',
                     loc=[-4.5020, -4.4628, -4.4237, -4.3846, -4.3455, -4.3064, -4.2672, -4.2281, -4.1890, -4.1499,
                          -4.1108, -4.0717, -4.0327, -3.9938, -3.9547, -3.9157, -3.8766, -3.8380, -3.7395, -3.6386,
                          -3.5423, -3.4537, -3.3799, -3.3195, -3.2744], locWidth=0.02):
    # loc = [-4.5020,   -4.4628,   -4.4237,   -4.3846,   -4.3455,   -4.3064,   -4.2672,   -4.2281, -4.1890,   -4.1499,   -4.1108,   -4.0717,   -4.0327,   -3.9938   -3.9547,   -3.9157, -3.8766,   -3.8380,   -3.7395,   -3.6386,   -3.5423,   -3.4537,   -3.3799,   -3.3195,-3.2744]
    # loc = [-4.50,-4.437,-4.3570,-4.2770,-4.1970,-4.1170,-4.0371,-3.9574,-3.8775,-3.5341,-3.3679,-3.273]
    # x1,x2,z1,z2,length,Z = plot2dGeom('input/iterRefinedTest.cfg')
    x1, x2, x3, y1, y2, y3, z1, z2, z3, area, Z, surfIndArray, surf, a, b, c, d, plane_norm = read3dGeom(
        'input/iterRefinedTest.cfg')
    plt.close()
    grossDep, grossEro, sumWeightStrike, E, A, EAdist, surfaceNumbers, surfEdist = nc_readSurface()
    plt.close()
    print('e', E)
    # loc1 = -4.17
    # loc2 = -4.081
    # loc3 = -4.25
    # loc4 = -3.6
    # locWidth = 0.02
    surfInd = surf > 0
    z1 = np.extract(surfInd, z1)
    gitrDirHe, heFlux = printHeDist(z=loc)
    os.mkdir('GITRoutput_W')
    nLocations = len(loc)
    for i in range(len(loc)):
        condition = [(z1 < loc[i] + locWidth) & (z1 > loc[i] - locWidth)]
    theseInd = np.where(condition)[1]
    print('theseInd', theseInd)
    print('condition', condition)
    print('surfEdist Shape', surfEdist.shape)
    print('theseInd', theseInd)
    print('theseInd', theseInd.shape)
    thisDist = surfEdist[theseInd, :, :]
    print('size thisDist', thisDist.shape)
    thisDist = np.sum(thisDist, axis=0)
    print('size thisDist', thisDist.shape)
    # thisDist = np.sum(thisDist,axis=1)
    print('size thisDist', thisDist.shape)
    print('thisDist', thisDist)
    dep = np.extract(condition, grossDep)
    ero = np.extract(condition, grossEro)
    strike = np.extract(condition, sumWeightStrike)
    areas = np.extract(condition, area)
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)
    # backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
    # backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
    # time = float(config.postProcessing.time);
    nParticles = float(config.impurityParticleSource.nP);
    # backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    erodedFlux = float(config.postProcessing.totalWFlux);
    erodedFluxPerParticle = erodedFlux / nParticles;
    netErosion = np.sum(ero - dep);
    netStrike = np.sum(strike)
    totalArea = np.sum(areas)
    impurityFlux = netErosion * erodedFluxPerParticle;
    if (heFlux[i] == 0.0):
        heFlux[i] = 1.0e19
    Wfrac = impurityFlux / heFlux[i];
    Aweight = np.sum(thisDist, axis=1)
    print('W impurity flux ', impurityFlux)
    print('W impurity fraction ', Wfrac)
    # for i in surfIndArray:
    gitrDir = 'GITRoutput_W/' + 'gitr' + str(i)
    os.mkdir(gitrDir)
    np.savetxt(gitrDir + '/gitrFluxE.dat', E)
    np.savetxt(gitrDir + '/gitrFluxAweight.dat', Aweight)
    np.savetxt(gitrDir + '/gitrFluxA.dat', A[:-1])
    np.savetxt(gitrDir + '/gitrFluxEAdist.dat', thisDist)

    if (path != ''):
        np.savetxt(path + '/' + gitrDir + '/gitrFluxE.dat', E)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxAweight.dat', Aweight)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxA.dat', A)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxEAdist.dat', thisDist)

    Dfrac = float(config.postProcessing.Dfrac);
    Hefrac = float(config.postProcessing.Hefrac);
    Tfrac = float(config.postProcessing.Tfrac);
    file = open('gitrOut_' + str(i) + '.txt', 'w')
    file.write('plasmaSpecies=He W D T\n')
    file.write('inputEnergy=-1.0 -1.0 0.0 0.0\n')
    file.write('inputAngle=-1.0 -1.0 0.0 0.0\n')
    file.write('fluxFraction=' + str(Hefrac) + ' ' + str(Wfrac) + ' ' + str(Dfrac) + ' ' + str(Tfrac) + '\n')
    file.write('flux=' + str(heFlux[i] / 1e18) + '\n')
    file.write('gitrOutputDir_He=' + os.getcwd() + '/' + 'GITRoutput/' + 'gitr' + str(i) + '\n')
    file.write('gitrOutputDir_W=' + os.getcwd() + '/' + gitrDir + '\n')
    file.close()

    if (path != ''):
        shutil.copyfile('gitrOut.txt', path + '/' + gitrDir + '/gitrOut.txt')
    # file = open(path+'/'+'gitrOut.txt','w')
    #
    # file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n')
    # file.write('flux='+str(backgroundFlux/1e18)+'\n')
    # file.write('gitrOutputDir='+os.getcwd()+'\n')
    # file.close()

    E0 = 0.0;
    E1 = 1000.0;
    nE = 200;
    dE = E1 / nE;
    Evec = np.linspace(0.5 * dE, E1 - 0.5 * dE, nE);
    A0 = 0.0;
    A1 = 90.0;
    nA = 30;
    for j in range(0, nA):
        # print('Evec size ', Evec.size)
        # print('EAdist[:,i] size ', EAdist[:,i].size)
        edOut = np.column_stack((Evec, EAdist[:, j]))
        np.savetxt(gitrDir + '/dist' + str(j) + '.dat', edOut)
    return nLocations


def iter3dProcessingQ4(nParticles=1200000, totalParticleRate=1.0, path='', locRmRs=[-0.1, 0.02, 0.09, 0.2],
                       flux_fracs=np.zeros((4, 5)), locWidth=0.02):
    x1, x2, x3, y1, y2, y3, z1, z2, z3, area, Z, surfIndArray, surf, a, b, c, d, plane_norm = read3dGeom(
        'input/iterRefinedTest.cfg')
    plt.close()
    grossDep, grossEro, sumWeightStrike, E, A, EAdist, surfaceNumbers, surfEdist = nc_readSurface()
    SOLPS = np.loadtxt('solpsTarg.txt', dtype='float', skiprows=1, delimiter=' ')
    rmrsSOLPS = SOLPS[:, 0]
    zSOLPS = SOLPS[:, 2]
    fluxInds = [28, 30, 32, 33, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
    fi = np.absolute(SOLPS[:, fluxInds]);
    totalFluxRmRs = np.sum(fi, axis=1)
    plt.close()
    print('e', E)
    surfInd = surf > 0
    z1 = np.extract(surfInd, z1)
    if os.path.isdir("GITRoutput_W"):
        print('GITRoutput_W already exists')
    else:
        os.mkdir('GITRoutput_W')
    frz = scii.interp1d(rmrsSOLPS, zSOLPS, fill_value=(zSOLPS[0], zSOLPS[-1]), bounds_error=False)
    fluxRmRsInterp = scii.interp1d(rmrsSOLPS, totalFluxRmRs, fill_value=(totalFluxRmRs[0], totalFluxRmRs[-1]),
                                   bounds_error=False)
    fluxes = fluxRmRsInterp(locRmRs)
    loc = frz(locRmRs)
    nLocations = len(loc)
    for i in range(len(loc)):
        condition = [(z1 < loc[i] + locWidth) & (z1 > loc[i] - locWidth)]
    theseInd = np.where(condition)[1]
    print('theseInd', theseInd)
    print('condition', condition)
    print('surfEdist Shape', surfEdist.shape)
    print('theseInd', theseInd)
    print('theseInd', theseInd.shape)
    thisDist = surfEdist[theseInd, :, :]
    print('size thisDist', thisDist.shape)
    thisDist = np.sum(thisDist, axis=0)
    print('size thisDist', thisDist.shape)
    # thisDist = np.sum(thisDist,axis=1)
    print('size thisDist', thisDist.shape)
    print('thisDist', thisDist)
    dep = np.extract(condition, grossDep)
    ero = np.extract(condition, grossEro)
    strike = np.extract(condition, sumWeightStrike)
    surfAreas = np.extract(surfInd, area)
    areas = np.extract(condition, surfAreas)
    print("number of areas", len(areas))
    print("areas", areas)

    # nParticles = float(config.impurityParticleSource.nP);
    # erodedFlux = float(config.postProcessing.totalWFlux);
    erodedFlux = totalParticleRate
    erodedFluxPerParticle = erodedFlux / nParticles;
    netErosion = np.mean(np.divide((ero - dep), areas));
    netStrike = np.sum(strike)
    totalArea = np.sum(areas)
    impurityFlux = erodedFluxPerParticle * np.mean(np.divide(strike, areas));
    print("loc netStrike totalArea", loc[i], netStrike, totalArea)
    # heatFluxRmRs = [-0.109154550719626 , -9.653786554458987E-002 , -7.327382054436374E-002 , -5.340755520166962E-002 , -3.684860232274187E-002 , -2.359568221752333E-002 , -1.342808911941418E-002 , -6.074684322167428E-003 , -1.331154862563744E-003 , 1.331154862563730E-003 , 4.767247358658980E-003 , 1.153375963886041E-002 , 2.146828469948289E-002 , 3.446721044754836E-002 , 5.053455048502209E-002 , 6.969946162664536E-002 , 9.215179987578859E-002 , 0.118084582163366 , 0.147824629366545 , 0.181757787874304 , 0.220748070819577 , 0.266226319146929 , 0.314553760432419 , 0.362245819191070 , 0.427758324849402 , 0.544486523268327 , 0.688286110504397 , 0.830542303703098 , 0.957948318740153 , 1.05083628019670 , 1.11603260806660 , 1.16594128670273 , 1.20417900393696 , 1.23354012318983 , 1.26056203356870 , 1.28515695577904 , 1.30459260804867 , 1.31243123811648]
    # heatFlux=[1.63487, 5377.34, 18518.8, 57355.4, 153011, 341486, 636782, 1.02231e+06, 1.4209e+06, 1.72315e+06, 2.17449e+06, 3.00281e+06, 3.84984e+06, 4.37564e+06, 5.56386e+06, 5.03941e+06, 3.6033e+06, 2.3909e+06, 1.59189e+06, 1.01545e+06, 595130, 335694, 226479, 144682, 52446.4, 18123.4, 12261.8, 6953.63, 5460.5, 4601.23, 3979.88, 2822.06, 2668.14, 1798.27, 1218.84, 750.364, 449.698, 4.43623]
    heatFluxRmRs = [-0.0968, -0.0735, -0.0537, -0.0371, -0.0239, -0.0138, -0.0064, -0.0017, 0.0009, 0.0044, 0.0112,
                    0.0211, 0.0341, 0.0502, 0.0693, 0.0918, 0.1177, 0.1475, 0.1814, 0.2204, 0.2658, 0.3142, 0.3619,
                    0.4274, 0.5442, 0.6889, 0.8319, 0.9599, 1.0529, 1.1182, 1.1682, 1.2065, 1.2358, 1.2629, 1.2874,
                    1.3069]

    heatFlux = np.multiply(1.0e+06,
                           [0.1106, 0.1359, 0.2066, 0.3500, 0.6497, 1.0946, 1.6078, 2.0966, 2.8286, 3.4223, 4.8043,
                            5.2989, 4.7962, 5.1764, 5.0547, 4.2162, 2.8670, 1.9490, 1.4311, 0.9886, 0.7624, 0.6522,
                            0.5521, 0.3994, 0.3265, 0.2717, 0.2252, 0.1756, 0.1473, 0.0932, 0.1185, 0.1028, 0.1000,
                            0.0798, 0.0781, 0.0142])
    hfInterp = scii.interp1d(heatFluxRmRs, heatFlux, fill_value=(heatFlux[0], heatFlux[-1]), bounds_error=False)
    thisHeatFlux = hfInterp(locRmRs[i])
    Wfrac = impurityFlux / fluxes[i];
    Aweight = np.sum(thisDist, axis=1)
    print('W impurity flux ', impurityFlux)
    print('W impurity fraction ', Wfrac)
    # for i in surfIndArray:
    gitrDir = 'GITRoutput_W/' + 'gitr' + str(i)
    if os.path.isdir(gitrDir):
        print(gitrDir + ' already exists')
    else:
        os.mkdir(gitrDir)
    np.savetxt(gitrDir + '/gitrFluxE.dat', E)
    np.savetxt(gitrDir + '/gitrFluxAweight.dat', Aweight)
    np.savetxt(gitrDir + '/gitrFluxA.dat', A[:-1])
    np.savetxt(gitrDir + '/gitrFluxEAdist.dat', thisDist)

    if (path != ''):
        np.savetxt(path + '/' + gitrDir + '/gitrFluxE.dat', E)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxAweight.dat', Aweight)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxA.dat', A)
        np.savetxt(path + '/' + gitrDir + '/gitrFluxEAdist.dat', thisDist)

    # Dfrac = float(config.postProcessing.Dfrac);
    # Hefrac = float(config.postProcessing.Hefrac);
    # Tfrac = float(config.postProcessing.Tfrac);
    file = open('gitrOut' + str(i) + '.txt', 'w')
    file.write('plasmaSpecies=He W D T Ne\n')
    file.write('fluxFraction=' + str(flux_fracs[0, i]) + ' ' + str(Wfrac) + ' ' + str(flux_fracs[1, i]) + ' ' + str(
        flux_fracs[2, i]) + ' ' + str(flux_fracs[3, i]) + '\n')
    file.write('flux=' + str(fluxes[i] / 1e18) + '\n')
    file.write('heat=' + str(thisHeatFlux / 1e18) + ' 343.0' + '\n')
    # file.write('gitrOutputDir_W='+os.getcwd()+'gitrOutput_W'+'\n')
    # file.write('gitrOutputDir_D='+os.getcwd()+'gitrOutput_D'+'\n')
    # file.write('gitrOutputDir_T='+os.getcwd()+'gitrOutput_T'+'\n')
    # file.write('gitrOutputDir_He='+os.getcwd()+'gitrOutput_He'+'\n')
    # file.write('gitrOutputDir_Ne='+os.getcwd()+'gitrOutput_Ne'+'\n')
    file.write('gitrOutputDir_W=' + os.getcwd() + '/GITRoutput_W/gitr' + str(i) + '\n')
    file.write('gitrOutputDir_D=' + os.getcwd() + '/GITRoutput_D/gitr' + str(i) + '\n')
    file.write('gitrOutputDir_T=' + os.getcwd() + '/GITRoutput_T/gitr' + str(i) + '\n')
    file.write('gitrOutputDir_He=' + os.getcwd() + '/GITRoutput_He/gitr' + str(i) + '\n')
    file.write('gitrOutputDir_Ne=' + os.getcwd() + '/GITRoutput_Ne/gitr' + str(i) + '\n')
    file.close()
    # file.write('inputEnergy=-1.0 -1.0 0.0 0.0\n')
    # file.write('inputAngle=-1.0 -1.0 0.0 0.0\n')
    # file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n')
    # file.write('flux='+str(heFlux[i]/1e18)+'\n')
    # file.write('gitrOutputDir_He='+os.getcwd()+'/'+'GITRoutput/'+'gitr'+str(i)+'\n')
    # file.write('gitrOutputDir_W='+os.getcwd()+'/'+gitrDir+'\n')
    # file.close()

    if (path != ''):
        shutil.copyfile('gitrOut.txt', path + '/' + gitrDir + '/gitrOut.txt')
    # file = open(path+'/'+'gitrOut.txt','w')
    #
    # file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n')
    # file.write('flux='+str(backgroundFlux/1e18)+'\n')
    # file.write('gitrOutputDir='+os.getcwd()+'\n')
    # file.close()

    E0 = 0.0;
    E1 = 1000.0;
    nE = 200;
    dE = E1 / nE;
    Evec = np.linspace(0.5 * dE, E1 - 0.5 * dE, nE);
    A0 = 0.0;
    A1 = 90.0;
    nA = 30;
    for j in range(0, nA):
        # print('Evec size ', Evec.size)
        # print('EAdist[:,i] size ', EAdist[:,i].size)
        edOut = np.column_stack((Evec, EAdist[:, j]))
        np.savetxt(gitrDir + '/dist' + str(j) + '.dat', edOut)

    return nLocations


def piscesProcessing(r=0.01, path=''):
    x1, x2, x3, y1, y2, y3, z1, z2, z3, area, Z, surfIndArray = read3dGeom('input/gitrGeometryPisces1inch.cfg')
    r1 = np.sqrt(np.multiply(x1[surfIndArray], x1[surfIndArray]) + np.multiply(y1[surfIndArray], y1[surfIndArray]))
    r2 = np.sqrt(np.multiply(x2[surfIndArray], x2[surfIndArray]) + np.multiply(y2[surfIndArray], y2[surfIndArray]))
    r3 = np.sqrt(np.multiply(x3[surfIndArray], x3[surfIndArray]) + np.multiply(y3[surfIndArray], y3[surfIndArray]))
    grossDep, grossEro, sumWeightStrike, E, A, EAdist = nc_readSurface()
    condition = r1 < r
    dep = np.extract(condition, grossDep)
    ero = np.extract(condition, grossEro)
    strike = np.extract(condition, sumWeightStrike)
    areas = np.extract(condition, area)
    with io.open('input/gitrInput.cfg') as f:
        config = libconf.load(f)

    backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec);  # 3.8640e+19;for pisces He high flux case
    backgroundFlux = float(config.postProcessing.backgroundFlux);  # 3.5e22;
    time = float(config.postProcessing.time);
    nParticles = float(config.impurityParticleSource.nP);
    backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
    erodedMass = time * backgroundIonsPerSec * 184 * 1.66e-27 * backgroundSputtYield * 1000;
    erodedMassPerParticle = erodedMass / nParticles;
    netErosion = np.sum(ero - dep);
    netStrike = np.sum(strike)
    totalArea = np.sum(areas)
    impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec * backgroundSputtYield / nParticles;
    impurityFlux = netStrike / totalArea * impurityParticlePerSecondPerComputationalPartice;
    Wfrac = impurityFlux / backgroundFlux;
    Aweight = np.sum(EAdist, axis=0)
    print('W impurity flux ', impurityFlux)
    print('W impurity fraction ', Wfrac)
    # for i in surfIndArray:
    np.savetxt('gitrFluxE.dat', E)
    np.savetxt('gitrFluxAweight.dat', Aweight)
    np.savetxt('gitrFluxA.dat', A[:-1])
    np.savetxt('gitrFluxEAdist.dat', EAdist)

    if (path != ''):
        np.savetxt(path + '/' + 'gitrFluxE.dat', E)
        np.savetxt(path + '/' + 'gitrFluxAweight.dat', Aweight)
        np.savetxt(path + '/' + 'gitrFluxA.dat', A)
        np.savetxt(path + '/' + 'gitrFluxEAdist.dat', EAdist)

    Dfrac = float(config.postProcessing.Dfrac);
    Hefrac = float(config.postProcessing.Hefrac);
    Tfrac = float(config.postProcessing.Tfrac);
    file = open('gitrOut.txt', 'w')
    file.write('plasmaSpecies=He W D T\n')
    file.write('fluxFraction=' + str(Hefrac) + ' ' + str(Wfrac) + ' ' + str(Dfrac) + ' ' + str(Tfrac) + '\n')
    file.write('flux=' + str(backgroundFlux / 1e18) + '\n')
    file.write('gitrOutputDir=' + os.getcwd() + '\n')
    file.close()

    if (path != ''):
        shutil.copyfile('gitrOut.txt', path + '/' + 'gitrOut.txt')
    # file = open(path+'/'+'gitrOut.txt','w')
    #
    # file.write('fluxFraction='+str(Hefrac)+' '+str(Wfrac)+ ' '+str(Dfrac)+ ' ' + str(Tfrac)+'\n')
    # file.write('flux='+str(backgroundFlux/1e18)+'\n')
    # file.write('gitrOutputDir='+os.getcwd()+'\n')
    # file.close()

    E0 = 0.0;
    E = 1000.0;
    nE = 100;
    dE = E / nE;
    Evec = np.linspace(0.5 * dE, E - 0.5 * dE, nE);
    A0 = 0.0;
    A = 90.0;
    nA = 90;
    for i in range(0, nA):
        print('Evec size ', Evec.size)
        print('EAdist[:,i] size ', EAdist[:, i].size)
        edOut = np.column_stack((Evec, EAdist[:, i]))
        np.savetxt('dist' + str(i) + '.dat', edOut)


def d3dProcessing(r=0.01, path=''):
    # plt.figure(1,figsize=(6, 10), dpi=1000)
    # plot2dGeom(filename='../input/gitrD3DGeometry2DWrings.cfg')
    grossDep, grossEro, sumWeightStrike, E, A, EAdist = nc_readSurface()
    print('W gross deposition ', grossDep)
    print('W gross erosion ', grossEro)


# condition = r1 < r
# dep = np.extract(condition,grossDep)
# ero = np.extract(condition,grossEro)
# strike = np.extract(condition,sumWeightStrike)
# areas = np.extract(condition,area)
# with io.open('input/gitrInput.cfg') as f:
#    config = libconf.load(f)

# backgroundIonsPerSec = float(config.postProcessing.backgroundIonsPerSec); #3.8640e+19;for pisces He high flux case
# backgroundFlux = float(config.postProcessing.backgroundFlux);#3.5e22;
# time = float(config.postProcessing.time);
# nParticles = float(config.impurityParticleSource.nP);
# backgroundSputtYield = float(config.postProcessing.backgroundSputtYield);
# erodedMass = time*backgroundIonsPerSec*184*1.66e-27*backgroundSputtYield*1000;
# erodedMassPerParticle = erodedMass/nParticles;
# netErosion = np.sum(ero - dep);
# netStrike = np.sum(strike)
# totalArea = np.sum(areas)
# impurityParticlePerSecondPerComputationalPartice = backgroundIonsPerSec*backgroundSputtYield/nParticles;
# impurityFlux = netStrike/totalArea*impurityParticlePerSecondPerComputationalPartice;
# Wfrac = impurityFlux/backgroundFlux;
# Aweight = np.sum(EAdist,axis=0)
# print('W impurity flux ', impurityFlux)
# print('W impurity fraction ', Wfrac)
def plot3dGeom(filename="gitrGeometry.cfg"):
    with io.open(filename) as f:
        config = libconf.load(f)

    x1 = np.array(config.geom.x1)
    x2 = np.array(config.geom.x2)
    x3 = np.array(config.geom.x3)
    y1 = np.array(config.geom.y1)
    y2 = np.array(config.geom.y2)
    y3 = np.array(config.geom.y3)
    z1 = np.array(config.geom.z1)
    z2 = np.array(config.geom.z2)
    z3 = np.array(config.geom.z3)
    xs = []
    ys = []
    zs = []
    for i in range(0, x1.size - 1):
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


# fig.patch.set_facecolor('black')
# ax = fig.gca(projection='3d')
# ax.plot_trisurf(xs,ys)
# print('xs ys zs', xs, ys, zs)
# ax = Axes3D(fig)
# verts = [zip(xs,ys,zs)]
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_xlim3d(0.5,2.5)
# ax.set_ylim3d(-0.2,0.2)
# ax.set_zlim3d(-0.5,1.5)
# ax.add_collection3d(Poly3DCollection(verts))
# plt.savefig('geomPlot.png')
def nc_plotHist(filename='output/history.nc', plot=1):
    ncFile = netCDF4.Dataset(filename, "r")
    nT = ncFile.dimensions['nT'].size
    nP = ncFile.dimensions['nP'].size
    x = np.reshape(np.array(ncFile.variables['x']), (nP, nT))
    y = np.reshape(np.array(ncFile.variables['y']), (nP, nT))
    z = np.reshape(np.array(ncFile.variables['z']), (nP, nT))
    vx = np.reshape(np.array(ncFile.variables['vx']), (nP, nT))
    vy = np.reshape(np.array(ncFile.variables['vy']), (nP, nT))
    vz = np.reshape(np.array(ncFile.variables['vz']), (nP, nT))
    charge = np.reshape(np.array(ncFile.variables['charge']), (nP, nT))
    weight = np.reshape(np.array(ncFile.variables['weight']), (nP, nT))
    r = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
    # print('z ', z[:,0])
    # print('x ', x.shape)
    # print('x ', x[0,:])
    # print('r ', r.shape)
    # print('z ', z[0][:].shape)
    ##for i in range(nT):
    ##    print('r,z ',i, x[i][0],z[i][0])
    # single = x[0][:];
    if plot == 1:
        plt.close()
        plt.figure(1, figsize=(6, 10), dpi=60)
        # plot2dGeom(filename='../input/gitrGeometry.cfg')
        if (x.shape[0] == 1):
            plt.plot(r[0][:], z[0][:], linewidth=5, color='green')
        else:
            for i in range(nP):
                # print('i', i)
                #              print('size', r[:,i].size)
                #              print('r ', r[:,i])
                plt.plot(r[i, :], z[i, :], linewidth=0.5)
        #              #plt.plot(r[i,:],z[i,:],linewidth=1.0)
        #              #plt.setp(linewidth=0.2)
        # plt.xlim((5.3,6.0))
        # plt.ylim((-4.4,-3.0))
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.title('DIII-D W Impurity Simulation', fontsize=20)
        plt.xlabel('r [m]', fontsize=16)
        plt.ylabel('z [m]', fontsize=16)
        # plt.axis('equal')
        print('saving tracksRZ')
        plt.savefig('tracksRZ.png')
        plt.show()
        plt.close()
        plt.figure(1, figsize=(10, 6), dpi=100)
        if (x.shape[0] == 1):
            plt.plot(x[0][:], y[0][:], linewidth=0.5)
        else:
            for i in range(nP):
                # print('i', i)
                plt.plot(x[i, :], y[i, :], linewidth=0.5)
        # plt.ylim((-0.02,0.02))
        # plt.xlim((-0.02,0.02))
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.axis('equal')
        print('saving tracksXY')
        plt.savefig('tracksXY.png')
    return x, y, z, r, vx, vy, vz, charge, weight, nP, nT


def nc_plotSpec(filename='spec.nc'):
    ncFile = netCDF4.Dataset(filename, "r")
    nBins = ncFile.dimensions['nBins'].size
    nR = ncFile.dimensions['nR'].size
    nZ = ncFile.dimensions['nZ'].size
    n = np.array(ncFile.variables['n'])
    print('n ', n.shape)
    dens = n[nBins - 1, :, :]
    print('dens ', dens.shape)
    plt.close()
    plt.figure(1, figsize=(10, 6), dpi=2000)
    plotsize = math.ceil(nBins ** (0.5))
    for i in range(nBins - 1, nBins):
        dens = np.log10(n[i, :, :])
        # plt.subplot(plotsize,plotsize,i+1)
        plot2dGeom('../2d/input/iter2dRefinedOuterTarget.cfg')
        plt.title("ITER W Impurity Density")
        plt.xlabel("r [m]")
        plt.ylabel("z [m]")
        plt.imshow(dens, extent=[4.0, 8.4, -4.6, 4.7], origin='lower')
        plt.colorbar(orientation='vertical')
    plt.savefig('image1.png')


def nc_plotSpec3D(filename='spec.nc'):
    ncFile = netCDF4.Dataset(filename, "r")
    nBins = ncFile.dimensions['nBins'].size
    nR = ncFile.dimensions['nR'].size
    nY = ncFile.dimensions['nY'].size
    nZ = ncFile.dimensions['nZ'].size
    n = np.array(ncFile.variables['n'])
    print('n ', n.shape)
    dens = n[nBins - 1, :, :, :]
    print('dens ', dens.shape)
    plt.close()
    plt.figure(1, figsize=(10, 6), dpi=2000)
    plt.imshow(n[0, 5, :, :])
    plt.colorbar(orientation='vertical')
    # plotsize = math.ceil(nBins**(0.5))
    # for i in range(nBins):
    #    dens = np.log10(n[i,:,:])
    #    plt.subplot(plotsize,plotsize,i+1)
    #    plt.imshow(dens,origin='lower')
    #    plt.colorbar(orientation='vertical')
    plt.savefig('slice1.png')


def nc_plotPositions(filename='output/positions.nc', cfg='input/gitrInput.cfg'):
    ncFile = netCDF4.Dataset(filename, "r")
    nP = ncFile.dimensions['nP'].size
    x = np.array(ncFile.variables['x'])
    y = np.array(ncFile.variables['y'])
    r = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
    z = np.array(ncFile.variables['z'])
    charge = np.array(ncFile.variables['charge'])
    print('x ', x)
    print('y ', y)
    print('r ', r)
    print('z ', z)
    plt.close()
    with io.open(cfg) as f:
        config = libconf.load(f)
    used3d = config['flags']['USE3DTETGEOM']
    hasTracks = config['flags']['PARTICLE_TRACKS']
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    f.suptitle('Position Plots for x,y(left) and r,z (right)')
    ax1.scatter(x, y)
    ax2.scatter(r, z)
    if (used3d == 0):
        gcfg = 'input/gitrGeometry.cfg'
        with io.open(gcfg) as ff:
            config = libconf.load(ff)

        x1 = np.array(config.geom.x1)
        x2 = np.array(config.geom.x2)
        z1 = np.array(config.geom.z1)
        z2 = np.array(config.geom.z2)
        ax2.plot(np.append(x1, x1[0]), np.append(z1, z1[0]), linewidth=2.0, color='k')
    if (hasTracks):
        ncFile = netCDF4.Dataset("output/history.nc", "r")
        nT = ncFile.dimensions['nT'].size
        nP = ncFile.dimensions['nP'].size
        xt = np.reshape(np.array(ncFile.variables['x']), (nP, nT))
        yt = np.reshape(np.array(ncFile.variables['y']), (nP, nT))
        zt = np.reshape(np.array(ncFile.variables['z']), (nP, nT))
        rt = np.sqrt(np.multiply(xt, xt) + np.multiply(yt, yt))
        plt.figure(1, figsize=(6, 10), dpi=60)
        if (xt.shape[0] == 1):
            ax1.plot(xt[0][:], yt[0][:], linewidth=1, color='green')
            ax2.plot(rt[0][:], zt[0][:], linewidth=1, color='green')
        else:
            for i in range(nP):
                ax1.plot(xt[i, :], yt[i, :], linewidth=0.5)
                ax2.plot(rt[i, :], zt[i, :], linewidth=0.5)

    # fig = plt.figure()
    # ax = fig.add_subplot(121)
    # ax.scatter(x, y)
    # ax = fig.add_subplot(122)
    # ax.scatter(r,z)
    # plt.title("Position Plots for x,y(left) and r,z (right)")
    ax1.set(xlabel='x-label', ylabel='y-label')
    ax2.set(xlabel='x-label', ylabel='y-label')
    # f.xlabel("r [m]")
    # f.ylabel("z [m]")
    f.savefig('positions.png')
    return x, y, r, z, charge


def nc_plotVz(filename='history.nc'):
    ncFile = netCDF4.Dataset(filename, "r")
    nT = ncFile.dimensions['nT'].size
    nP = ncFile.dimensions['nP'].size
    x = np.reshape(np.array(ncFile.variables['x']), (nT, nP))
    y = np.reshape(np.array(ncFile.variables['y']), (nT, nP))
    z = np.reshape(np.array(ncFile.variables['z']), (nT, nP))
    vx = np.reshape(np.array(ncFile.variables['vx']), (nT, nP))
    vy = np.reshape(np.array(ncFile.variables['vy']), (nT, nP))
    vz = np.reshape(np.array(ncFile.variables['vz']), (nT, nP))
    vperp = np.sqrt(np.multiply(vx[nT - 1, :], vx[nT - 1, :]) + np.multiply(vy[nT - 1, :], vy[nT - 1, :]))
    pitchAngle = np.arctan(np.divide(vperp, vz[nT - 1, :]))
    r = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
    plt.figure(1, figsize=(10, 6), dpi=250)
    # plt.plot(z[:,1],vz[:,1],linewidth=0.5)
    plt.plot(z, vz, linewidth=0.5, label=str(np.linspace(0, nP - 1, nP)))
    # plt.legend()
    # plt.savefig('vz.png')
    # for i in range(nP):
    # print('i', i)
    # if((i > 20) and (i < 41)):
    # plt.plot(z[:,i],vz[:,i],linewidth=0.5,label=str(i*90.0/180))
    # print('size', r[:,i].size)
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
    # plt.legend()
    plt.savefig('vz.png')
    plt.close()
    plt.figure(1, figsize=(10, 6), dpi=250)
    plt.hist(pitchAngle, bins=30)
    plt.savefig('pa.png')


def nc_readSurface(filename='output/surface.nc'):
    import os
    ncFile = netCDF4.Dataset(filename, "r")
    nS = len(ncFile.dimensions['nSurfaces'])
    nE = len(ncFile.dimensions['nEnergies'])
    nA = len(ncFile.dimensions['nAngles'])
    E = np.linspace(0, 1000.0, nE + 1)
    A = np.linspace(0, 90.0, nA + 1)
    surfaceNumbers = np.array(ncFile.variables['surfaceNumber'])
    grossDep = np.array(ncFile.variables['grossDeposition'])
    grossEro = np.array(ncFile.variables['grossErosion'])
    sumWeightStrike = np.array(ncFile.variables['sumWeightStrike'])
    surfEdist = np.reshape(np.array(ncFile.variables['surfEDist']), (nS, nA, nE))
    surfEdist.reshape((nS, nA, nE))
    print('size surfEdist ', surfEdist.size)
    print('shape surfEdist ', surfEdist.shape)
    EAdist = np.sum(surfEdist, axis=0)
    print('shape EAdist ', EAdist.shape)
    EAdist = EAdist.reshape(nE, nA)
    plt.pcolor(EAdist)
    plt.colorbar()
    plt.savefig('EAdist.png')
    return grossDep, grossEro, sumWeightStrike, E, A, EAdist, surfaceNumbers, surfEdist


def plotPitch(filename='positions.nc'):
    ncFile = netCDF4.Dataset(filename, "r")
    nP = ncFile.dimensions['nP'].size
    x = np.array(ncFile.variables['x'])
    y = np.array(ncFile.variables['y'])
    z = np.array(ncFile.variables['z'])
    vx = np.array(ncFile.variables['vx'])
    vy = np.array(ncFile.variables['vy'])
    vz = np.array(ncFile.variables['vz'])
    vperp = np.sqrt(np.multiply(vx, vx) + np.multiply(vy, vy))
    pitchAngle = np.arctan(np.divide(vperp, vz))
    r = np.sqrt(np.multiply(x, x) + np.multiply(y, y))
    plt.figure(1, figsize=(10, 6), dpi=250)
    plt.hist(pitchAngle, bins=30)
    plt.savefig('pa.png')
    plt.close()
    plt.subplot(3, 1, 1)
    plt.hist(vx, bins=30)
    plt.subplot(3, 1, 2)
    plt.hist(vy, bins=30)
    plt.subplot(3, 1, 3)
    plt.hidast(vz, bins=30)
    plt.savefig('vs.png')


def make_gitr_geometry_from_solps(gitr_geometry_filename='gitr_geometry.cfg', \
                                  solps_mesh_extra='/Users/tyounkin/Downloads/mesh.extra.iter', \
                                  right_target_file='/Users/tyounkin/Code/solps-iter-data/build/rightTargOutput', \
                                  left_target_file='/Users/tyounkin/Code/solps-iter-data/build/leftTargOutput'):
    # This program uses the solps-iter mesh.extra file in combination
    # with the inner and outer (left and right) divertor target
    # coordinates which come from the solps-iter-data interpolation
    # program to create a 2d geometry for GITR in which
    # the solps plasma profiles properly match the divertor target
    # geometry.
    #
    # This geometry is then written to a config (cfg) file for
    # use in GITR simulation

    solps_mesh = np.loadtxt(solps_mesh_extra)

    plt.plot(solps_mesh[:, [0, 2]].transpose(), solps_mesh[:, [1, 3]].transpose())
    plt.show()

    manual_indices = np.array(range(88, 98))
    manual_indices = np.append(manual_indices, range(71, 82))
    manual_indices = np.append(manual_indices, range(68, 71))
    manual_indices = np.append(manual_indices, 50)
    manual_indices = np.append(manual_indices, range(36, 45))
    manual_indices = np.append(manual_indices, [34, 52, 51, 35])
    manual_indices = np.append(manual_indices, range(0, 34))
    manual_indices = np.append(manual_indices, 49)
    manual_indices = np.append(manual_indices, range(53, 68))
    manual_indices = np.append(manual_indices, range(82, 87))
    manual_indices = np.append(manual_indices, [127, 87])
    # [89:98,72:82,69:71,51,37:45,35,53,52,36,1:34,50,54:68
    # ,83:87,128,88]; Matlab indices

    print(manual_indices)
    r_iter = solps_mesh[manual_indices, 0]
    z_iter = solps_mesh[manual_indices, 1]
    r_iter = np.append(r_iter, solps_mesh[87, 2])
    z_iter = np.append(z_iter, solps_mesh[87, 3])

    r_dome = solps_mesh[range(98, 117), 0]
    z_dome = solps_mesh[range(98, 117), 1]
    r_dome = np.append(r_dome, solps_mesh[134, 0])
    z_dome = np.append(z_dome, solps_mesh[134, 1])

    r_iter = np.append(r_iter, r_dome)
    z_iter = np.append(z_iter, z_dome)

    plt.plot(r_iter, z_iter)
    plt.scatter(r_iter,z_iter)
    print('manual geometry size',r_iter.size)

    right_target = np.loadtxt(right_target_file)
    r_right_target = right_target[:, 0]
    z_right_target = right_target[:, 1]
    print('right target size',r_right_target.size)
    plt.scatter(r_right_target,z_right_target)

    left_target = np.loadtxt(left_target_file)
    r_left_target = left_target[:, 0]
    z_left_target = left_target[:, 1]
    print('left target size',r_left_target.size)
    plt.scatter(r_left_target,z_left_target)
    plt.show()

    print('r_final size before',r_iter.size)
    r_final, z_final = replace_line_segment(r_right_target, z_right_target, r_iter, z_iter)
    print('r_final size after',r_final.size)
    r_final, z_final = replace_line_segment(r_left_target, z_left_target, r_final, z_final)
    print('r_final size after',r_final.size)
    print(r_final.size)

    Z = np.zeros(len(r_final)+1)
    surfaces = np.zeros(len(r_final)+1)
    inDir = np.zeros(len(r_final))

    i_a, i_b = intersection(r_final, z_final, r_right_target, z_right_target)

    Z[i_b] = 74;
    surfaces[i_b] = 1;

    i_a, i_b = intersection(r_final, z_final, r_left_target, z_left_target)
    print('i_a',i_a)
    Z[i_b] = 74;
    surfaces[i_b] = 1;
    inDir[i_b] = -1;


    lines = gitr_lines_from_points(r_final, z_final)

    lines_to_gitr_geometry(gitr_geometry_filename, lines, Z, surfaces, inDir)

    removeQuotes(infile=gitr_geometry_filename, outfile=gitr_geometry_filename+"0")

    remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)
def removeQuotes(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace('"', ''))

def remove_endline_after_comma(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',\n', ','))

def remove_endline_after_comma2(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',       ', ','))

def lines_to_gitr_geometry(filename, lines, Z, surface, inDir):

    x1 = lines[:, 0]
    z1 = lines[:, 1]
    x2 = lines[:, 2]
    z2 = lines[:, 3]
    slope = lines[:, 4]
    intercept = lines[:, 5]
    line_length = lines[:, 6]
    fileExists = os.path.exists(filename)

    if not fileExists:
        f = open(filename, "w")
        f.close()

    with io.open(filename) as f:
        config = libconf.load(f)

    config['geom'] = {}
    config['geom']['x1'] = x1.tolist()
    config['geom']['z1'] = z1.tolist()
    config['geom']['x2'] = x2.tolist()
    config['geom']['z2'] = z2.tolist()
    config['geom']['slope'] = ['%.6f' % elem for elem in slope.tolist()]
    config['geom']['intercept'] = ['%.6f' % elem for elem in intercept.tolist()]
    config['geom']['length'] = line_length.tolist()
    config['geom']['Z'] = Z.tolist()
    config['geom']['surface'] = ['%i' % elem for elem in surface.tolist()]
    config['geom']['inDir'] = ['%i' % elem for elem in inDir.tolist()]
    config['geom']['y1'] = 0.0
    config['geom']['y2'] = 0.0
    config['geom']['periodic'] = 0

    with io.open(filename, 'w') as f:
        libconf.dump(config, f)


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


def replace_line_segment(x_priority, y_priority, x_base, y_base):
    x_final = x_base;
    y_final = y_base;

    all_close = np.array(0);
    for i in range(len(x_priority) - 1):
        distances = np.sqrt(np.power(x_priority[i] - x_base, 2) + np.power(y_priority[i] - y_base, 2));
        condition = (distances < np.sqrt(
            np.power((x_priority[i] - x_priority[i + 1]), 2) + np.power((y_priority[i] - y_priority[i + 1]), 2)));
        close_ones = np.where(condition)[0];
        all_close = np.append(all_close,close_ones);

    all_close = np.delete(all_close,0)
    remove_indices = np.unique(all_close);
    print('remove indices',remove_indices)
    d1 = np.sqrt(
        np.power((x_priority[0] - x_base[remove_indices[0]]), 2) + np.power((y_priority[0] - y_base[remove_indices[0]]),
                                                                            2));
    d2 = np.sqrt(np.power((x_priority[0] - x_base[remove_indices[-1]]), 2) + np.power(
        (y_priority[0] - y_base[remove_indices[-1]]), 2));

    if (d2 < d1):
        x_priority = np.flip(x_priority, 0)
        y_priority = np.flip(y_priority, 0)

    print('x_priority',x_priority)

    x_final = np.append(x_base[0:(remove_indices[0] - 1)], x_priority)
    x_final = np.append(x_final, x_base[(remove_indices[-1] + 1):]);
    y_final = np.append(y_base[0:(remove_indices[0] - 1)], y_priority)
    y_final = np.append(y_final, y_base[(remove_indices[-1] + 1):]);

    return x_final, y_final

def gitr_lines_from_points(r,z):

    nPoints = len(r);
    lines = np.zeros([nPoints, 7]);
    lines[:, 0] = r;
    lines[:, 1] = z;
    lines[0:-2, 2] = r[1:-1];
    lines[0:-2, 3] = z[1: -1];
    lines[-1, 2] = r[0];
    lines[-1, 3] = z[0];

    tol = 1e12;
    tol_small = 1e-12;

    for i in range(nPoints):
        if (lines[i, 3] - lines[i, 1]) == 0:
            lines[i, 4] = 0;
            lines[i, 5] = lines[i, 1];
        elif ((lines[i, 2] - lines[i, 0]) == 0):
            lines[i, 4] = np.sign(lines[i, 3] - lines[i, 1]) * tol;
            lines[i, 5] = tol;

        else:
            lines[i, 4] = (lines[i, 3] - lines[i, 1]) / (lines[i, 2] - lines[i, 0]);
            lines[i, 5] = -lines[i, 4] * lines[i, 0] + lines[i, 1];


    lines[:, 6] = np.sqrt((lines[:, 2] - lines[:, 0])**2 + (lines[:, 3] - lines[:, 1])** 2);

    return lines
if __name__ == "__main__":
    make_gitr_geometry_from_solps()
# asdfanc_show("surface.nc")
# depositedEdist()
# if(os.path.exists('output/history.nc')):
# 	nc_plotHist('output/history.nc')
# if(os.path.exists('output/spec.nc')):
#	nc_plotSpec('output/spec.nc')
# iter2dProcessing()
# iter3dProcessingQ4()
# printHeDist()
# nc_plotSpec3D()
# nc_plotPositions()
# nc_plotVz()
# plotPitch()
# piscesProcessing()
# modifyInputParam()
# nc_readSurface()
