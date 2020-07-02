#Author: Alyssa Hayes, University of Tennessee

import netCDF4
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator as interp
import matplotlib.pyplot as plt
from math import pi, exp, sqrt, log

#set global constants
nP = int(1e6)
k = 1.602e-19
gamma = 7
ke0 = 2000
sinj = 0.15
sv = 1.2

#---------------------------------------------------------

#initiate netcdf file
def init_ncid(case, fileExtension, L):
    nR = 10;
    r = np.linspace(-300,300,nR)
    nZ = 101
    z = np.linspace(0,2*L,nZ)
    mid = int(nZ/2)

    ncid = netCDF4.Dataset('../examples/SFT/%s/input/profiles'%case+fileExtension+'.nc','w')
    dimR = ncid.createDimension('nR',nR)
    dimZ = ncid.createDimension('nZ',nZ)
    ncid.createVariable('gridR','d',('nR',))[:] = r
    ncid.createVariable('gridZ','d',('nZ',))[:] = z
    return nR, r, nZ, z, mid, ncid

#assume Ti profile can be estimated by modified heat conduction equation (Eq. 26 in Stangeby and Elder, 1995)
def init_Ti(nR, z, mid, TD0, gamma, n0, fcond, ke0):
    #assume Ti=Te and determine the most probable speed from the initial T
    cs0 = sqrt((2*k*TD0) / (2*1.66e-27))

    #parallel power flux density at target
    P = gamma*n0*cs0*k*TD0

    #spatial multiplication factor in Eq. 26 for simplicity
    a = (7/2)*fcond*P / (ke0*(TD0**(7/2)))

    #define Ti(s) from Eq. 26
    Ti = TD0 * ((1+a*z)**(2/7))
    Ti[mid:] = np.flip(Ti[:mid+1])
    Ti2D = np.tile(Ti, [nR,1])

    #set inital distribution for rate of temperature change
    gradTi = (2*TD0*a) / (7*((a*z+1)**(5/7)))
    gradTi[mid:] = -np.flip(gradTi[:mid+1])
    gradTi2D = np.tile(gradTi, [nR,1])

    return cs0, Ti2D, gradTi2D

#extrapolate background flow velocity to a step function profile
def init_v(L, sv, nR, z, nZ, vb):
    #update vz due to vb
    vz = np.zeros(nZ)
    for i in np.where(z-sv<=np.finfo(float).eps): vz[i] = -vb
    for i in np.where(z+sv-2*L>=(np.finfo(float).eps)): vz[i] = vb
    vz2D = np.tile(vz, [nR,1])
    return vz2D

#populate netcdf files
def pop_ncid(ncid, vz2D, Ti2D, gradTi2D):
    ncid.createVariable('vx','d',('nR','nZ',))[:] = 0*vz2D
    ncid.createVariable('vy','d',('nR','nZ',))[:] = 0*vz2D
    ncid.createVariable('vz','d',('nR','nZ',))[:] = vz2D
    ncid.createVariable('ti','d',('nR','nZ',))[:] = Ti2D
    ncid.createVariable('te','d',('nR','nZ',))[:] = Ti2D
    ncid.createVariable('gradTiR','d',('nR','nZ',))[:] = 0*gradTi2D
    ncid.createVariable('gradTeR','d',('nR','nZ',))[:] = 0*gradTi2D
    ncid.createVariable('gradTiY','d',('nR','nZ',))[:] = 0*gradTi2D
    ncid.createVariable('gradTeY','d',('nR','nZ',))[:] = 0*gradTi2D
    ncid.createVariable('gradTiZ','d',('nR','nZ',))[:] = gradTi2D
    ncid.createVariable('gradTeZ','d',('nR','nZ',))[:] = 0*gradTi2D
    ncid.close()

#---------------------------------------------------------

#define profiles
def get_profilesA(case='caseA', fileExtension='', vb=750):
    L=10
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #define constants
    TD0 = 10
    dTds = 1.38

    #get velocity profile
    vz2D = init_v(L, sv, nR, z, nZ, vb)

    #set initial temperature distribution
    Ti = np.zeros(nZ)
    Ti[:mid] = TD0 + dTds*z[:mid]
    Ti[mid:] = TD0 + dTds*z[mid] + dTds*(L-z[mid:])
    Ti2D = np.tile(Ti, [nR,1])

    #set initial distribution for rate of temperature change
    gradTi = np.zeros(nZ)
    gradTi[:mid+1] = dTds
    gradTi[mid:] = -dTds
    gradTi2D = np.tile(gradTi, [nR,1])

    #export to ncid file
    pop_ncid(ncid, vz2D, Ti2D, gradTi2D)

def get_profilesB(case='caseB', fileExtension='', TD0=100):
    L = 30
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #define constants
    n0 = 1e20
    fcond = 0.1

    #get velocity, temperature, and temperature gradient profiles
    cs0, Ti2D, gradTi2D = init_Ti(nR, z, mid, TD0, gamma, n0, fcond, ke0)
    vz2D = init_v(L, sv, nR, z, nZ, 0.1*cs0)

    #export to ncid file
    pop_ncid(ncid, vz2D, Ti2D, gradTi2D)

def get_profilesC(case='caseC', fileExtension=''):
    L = 10
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #update electric field
    Ez0 = -25
    Ez = np.zeros(nZ)
    for i in np.where(z<=sv): Ez[i] = Ez0
    Ez2D = np.tile(Ez, [nR,1])

    #populate netcdf file
    ncid.createVariable('Ex','d',('nR','nZ',))[:] = 0*Ez2D
    ncid.createVariable('Ey','d',('nR','nZ',))[:] = 0*Ez2D
    ncid.createVariable('Ez','d',('nR','nZ',))[:] = Ez2D
    ncid.close()

def get_profilesD(case='caseD', fileExtension=''):
    L = 10
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #update vz from vb
    vb = -112.2
    vz = np.zeros(nZ)
    for i in np.where(z-sv<=np.finfo(float).eps): vz[i] = vb
    vz2D = np.tile(vz, [nR,1]).transpose()

    #populate netcdf file
    ncid.createVariable('vx','d',('nZ','nR',))[:] = 0*vz2D
    ncid.createVariable('vy','d',('nZ','nR',))[:] = 0*vz2D
    ncid.createVariable('vz','d',('nZ','nR',))[:] = vz2D
    ncid.close()

def get_profilesE(case='caseE', fileExtension='', TD0=10):
    L = 10
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #define constants
    n0 = 1e19
    fcond = 1

    #get velocity, temperature, and temperature gradient profiles
    cs0, Ti2D, gradTi2D = init_Ti(nR, z, mid, TD0, gamma, n0, fcond, ke0)
    vz2D = init_v(L, sv, nR, z, nZ, cs0)

    #export to ncid file
    pop_ncid(ncid, vz2D, Ti2D, gradTi2D)

def get_profilesF(case='caseF', fileExtension='', TD0=20):
    L = 30
    sv = 10.85
    nR, r, nZ, z, mid, ncid = init_ncid(case, fileExtension, L)

    #define constants
    n0 = 1e19
    fcond = 1

    #get velocity, temperature, and temperature gradient profiles
    cs0, Ti2D, gradTi2D = init_Ti(nR, z, mid, TD0, gamma, n0, fcond, ke0)
    vz2D = init_v(L, sv, nR, z, nZ, 0.1*cs0)

    #export to ncid file
    pop_ncid(ncid, vz2D, Ti2D, gradTi2D)

#---------------------------------------------------------

#generate particle source
def get_particleSource(fileExtension='', T=10, L=10):
    #set bounds
    z0 = 0.1
    z1 = 0.2
    x0 = -0.02
    x1 = 0.02

    #set mirrored particle source locations
    x = (x1-x0) * np.random.rand(nP) + x0
    y = np.zeros(nP)
    z = (z1-z0) * np.random.rand(nP) + z0
    odds = np.arange(1,nP,2)
    for i in np.nditer(odds): z[i] = z[i] - 2*(z1-z0) + 2*L - z0

    #set particle source velocity distribution
    k = 1.38e-23*11604
    B = 1.66e-27*12 / (2*k*T)
    vth = sqrt(2*k*T / (1.66e-27*12))
    vgrid = np.linspace(-3*vth, 3*vth)
    fv1 = sqrt(B/pi) * np.exp(-B*np.square(vgrid))
    fv1CDF = np.cumsum(fv1)/(np.cumsum(fv1)[-1])

    #interpolate from velocity distribution
    vx = scipy.interpolate.PchipInterpolator(fv1CDF, vgrid)(np.random.rand(nP))
    vy = scipy.interpolate.PchipInterpolator(fv1CDF, vgrid)(np.random.rand(nP))
    vz = scipy.interpolate.PchipInterpolator(fv1CDF, vgrid)(np.random.rand(nP))

    #create netcdf file
    ncid = netCDF4.Dataset('input/particleSource'+fileExtension+'.nc','w')
    dimP = ncid.createDimension('nP',nP)
    ncid.createVariable('x','d',('nP',))[:] = x
    ncid.createVariable('y','d',('nP',))[:] = y
    ncid.createVariable('z','d',('nP',))[:] = z
    ncid.createVariable('vx','d',('nP',))[:] = vx
    ncid.createVariable('vy','d',('nP',))[:] = vy
    ncid.createVariable('vz','d',('nP',))[:] = vz
    ncid.close()




#get characteristic slowing-down time
def drag_tau(mi, Zi, mD, TD, nD):
    nD = nD/1e18
    lnGam = 15
    top = mi*TD*sqrt(TD/mD)
    bottom = 6.8e4 * (1+mD/mi) * nD*lnGam*Zi**2
    return top/bottom



def SFTD(n0, gridZ, ne=1e19, TD=10, vD=-112.2):
    #define constants
    Z = 4
    m = 12
    tau_s = drag_tau(12,4,2,TD,ne)
    vth = sqrt(1.602e-19*TD/(12*1.66e-27))
    D_par = tau_s*vth**2
    vdiff = D_par/sinj
    vpl = vdiff + vD
    phi_in = 1.73e23
    nP = phi_in/vpl

    #get frictional force
    s = np.linspace(0,6,100)
    FFf = (1.66e-27*m*vD/(1.602e-19*TD*tau_s)) * gridZ

    #get impurity density
    ns = n0*np.exp(FFf)
    return ns



#compare impurity density decay rate to SFT
def SFTanalysis(cases=['D']):
    for i in range(0,len(cases)):

        #get outputs

        #generate output data
        outputDir='/Users/Alyssa/Dev/GITR/examples/SFT/case%s/output' %cases[i]
        specFilepath = '%s/spec1.nc' %outputDir
        specFile = netCDF4.Dataset(specFilepath,'r')
        nZ = specFile.dimensions['nZ'].size
        gridZ = np.array(specFile.variables['gridZ'])[:int(nZ/2)]
        dens = np.array(specFile.variables['n'])

        #process outputs to get impurity density profile
        c4dens = dens[4,:,:]
        nsum = 1.6586e12*np.sum(c4dens,1)[:int(nZ/2)]

        #reduce data to density decay region
        nmax = np.argmax(nsum)
        nreduced = np.trim_zeros(nsum)[nmax:]
        gridZred = gridZ[nmax:nmax+len(nreduced)]
        nsumred = nsum[nmax:nmax+len(nreduced)]

        #analyze data

        #fit exponential trend to decay rate
        trend_nsum = np.polyfit(gridZred, np.log(nreduced), 1)
        nfit = np.poly1d(trend_nsum)
        nfit0 = np.exp(nfit(gridZred)[0])

        #mean square error check
        diff = (nsumred - np.exp(nfit(gridZred)))/1e10
        rmse = np.sqrt(np.sum(np.square(diff))/len(diff))*1e10
        norm_factor = np.sum(np.exp(nfit(gridZred)))/len(gridZred)
        print('Root Mean Square Error = ', rmse/norm_factor)

        #get SFT results
        SFTns = SFTD(nfit0, gridZred)

        #get trendline of raw data ratio (slope of line in 2nd plot)
        trend_nsum_v_SFT = np.polyfit(SFTns, nsumred, 1)
        nsum_v_SFT_fit = np.poly1d(trend_nsum_v_SFT)
        GvS_slope = nsum_v_SFT_fit.coeffs[1]
        norm_GvS_slope = GvS_slope / (np.sum(SFTns)/len(SFTns))
        print('Normalized Density Ratio = ', norm_GvS_slope)

        #generate plots

        #plot over gridZ
        plt.semilogy(gridZ, nsum, gridZred, np.exp(nfit(gridZred)),gridZred,SFTns)
        plt.title('C4+ Impurity Density Profile')
        plt.xlabel('s[m]')
        plt.ylabel('n [m^-3]')
        plt.legend(('GITR Data','GITR Trendline','SFT Prediction'))
        plt.show()

        #plot GITR predictions over SFT predictions
        plt.plot(SFTns, nsumred,'.', SFTns, nsum_v_SFT_fit(SFTns))
        plt.title('GITR and SFT Density Comparison')
        plt.xlabel('SFT C4+ Density')
        plt.ylabel('GITR C4+ Density')
        plt.show()

if __name__ == "__main__":
    get_profilesF()
