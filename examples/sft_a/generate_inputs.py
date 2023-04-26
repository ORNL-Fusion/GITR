import numpy as np
import numpy.matlib
import netCDF4
#generate particle source and profiles

##################
# INPUTS SECTION #
##################

# Simulated particle mass (amu)
m=12;

# Domain half-length
L=10; # meters

# Flow speed toward the left boundary
flowSpeed = -750;

# Location at which plasma flow ends
sv = 1.2; #meters

# Ion temperature at the left boundary
TD0 = 10;

# Temperature gradient
dT_ds = 1.38;

# Number of simulated particles
nP=1e6;

# Spatial particle distribution - evenly between 10 and  20 cm
z0=0.1;
z1=0.2;
x0 = 0;
x1= 0;

# Plasma profile resolution
nR = 10;
nZ = 10001;

#########################
# END OF INPUTS SECTION #
#########################

##################################
# DEFINITIION OF PARTICLE SOURCE #
##################################

# Initialize particle positions randomly between x0 to x1, z0 to z1, y=0
x = (x1-x0)*np.random.uniform(0,1,int(nP))+x0;
z = (z1-z0)*np.random.uniform(0,1,int(nP))+z0;
y = 0.0*np.ones(int(nP));

# Initialize particle velocities to be thermal distribution at TD0
# Create 1D Maxwellian spped distribution
vTh = np.sqrt(2*TD0*1.602e-19/m/1.66e-27);
k = 1.38e-23*11604;
B = m*1.66e-27/(2*TD0*k);
vgrid = np.linspace(-3*vTh,3*vTh,1000000);
fv1 = np.sqrt(B/np.pi)*np.exp(-B*np.multiply(vgrid,vgrid));
fv1CDF = np.cumsum(fv1);
fv1CDF = fv1CDF/fv1CDF[-1];

# Sample 1D Maxwellian speed distribution in all velocity components
vx = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);
vy = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);
vz = np.interp(np.random.uniform(0,1,int(nP)),fv1CDF,vgrid);

# Write to netcdf file
rootgrp = netCDF4.Dataset("input/particleSource.nc", "w", format="NETCDF4")
npp = rootgrp.createDimension("nP", int(nP))
xxx = rootgrp.createVariable("x","f8",("nP"))
yyy = rootgrp.createVariable("y","f8",("nP"))
zzz = rootgrp.createVariable("z","f8",("nP"))
vxx = rootgrp.createVariable("vx","f8",("nP"))
vyy = rootgrp.createVariable("vy","f8",("nP"))
vzz = rootgrp.createVariable("vz","f8",("nP"))
xxx[:] = x
yyy[:] = y
zzz[:] = z
vxx[:] = vx
vyy[:] = vy
vzz[:] = vz
rootgrp.close()

##########################
# END OF PARTICLE SOURCE #
##########################

##################################
# DEFINITIION OF PLASMA PROFILES #
##################################

# Define 2D plasma profile grid
r = np.linspace(-300,300,nR);
z = np.linspace(0.0,2*L,nZ);

# Flow velocity setup
vz = np.zeros(nZ);
vz[np.where(z<=sv)] = flowSpeed;
vz2D = np.matlib.repmat(vz,nR,1);

# Ion temperature
Ti =np.zeros(nZ);
Ti = TD0+dT_ds*z;
Ti2D = np.matlib.repmat(Ti,nR,1);

vz2D = np.transpose(vz2D)
Ti2D = np.transpose(Ti2D)

# Write 2D plasma profiles to netcdf
rootgrp = netCDF4.Dataset("input/profiles.nc", "w", format="NETCDF4")
nrr = rootgrp.createDimension("nR", int(nR))
nzz = rootgrp.createDimension("nZ", int(nZ))

rrr = rootgrp.createVariable("gridR","f8",("nR"))
zzz = rootgrp.createVariable("gridZ","f8",("nZ"))
vzz = rootgrp.createVariable("vz","f8",("nZ","nR"))
vyy = rootgrp.createVariable("vy","f8",("nZ","nR"))
vxx = rootgrp.createVariable("vx","f8",("nZ","nR"))
tii = rootgrp.createVariable("ti","f8",("nZ","nR"))
tee = rootgrp.createVariable("te","f8",("nZ","nR"))
rrr[:] = r
zzz[:] = z
vxx[:] = 0*vz2D
vyy[:] = 0*vz2D
vzz[:] = vz2D
tii[:] = Ti2D
tee[:] = Ti2D
rootgrp.close()

##########################
# END OF PARTICLE SOURCE #
##########################
