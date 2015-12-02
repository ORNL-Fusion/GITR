%input file for GITR0.0
%xmin
xMinV =-.01;
% xmax
xMaxV =0.0;
% ymin
yMin = -0.02;
% ymax
yMax = 0.02;
% zmin
zMin = -0.02;
% zmax
zMax = 0.02;
% Surface cells in y
nY = 30;
% Surface cells in z
nZ = 40;
% Volume cells in x
nXv = 100;
% Volume cells in y
nYv = 80;
% Volume cells in z
nZv = 50;
% Sheath potential
sheathPotential = -60;
% Potential decay length
sheathWidth = 0.0001;
% Bfieldx
Bfieldx_dat =0.1;
% Bfieldy
Bfieldy_dat =1.0;
% Bfieldz
Bfieldz_dat =-0.3;
% Density
maxDensity = 1e19;
% Density Decay length
densitySOLDecayLength =0.1;
%Temperature (in eV)
maxTemp_eV = 20;
%Temperature decay length
tempSOLDecayLength = .1;
%Dperp
Dperp_dat = 0.04;
%Number of particles
nP = 1;
%Starting position in x
x_dat = -0.0075;
%Starting position in y
y_dat = -0.01;
%Starting position in z
z_dat = 0.0;
%Starting energy in x
Ex_dat = 1.0;
%Starting energy in y
Ey_dat = 1.5;
%Starting energy in z
Ez_dat = -1.0;
%Particle Mass
mass_dat = 12.0;
%Particles Charge State
charge_dat = 1.0;
%Slope of parameterized surface (X = m*Y +b)
surf_slope = -0.1;
%Intercept of parameterized surface
surf_incpt = -0.002;
%Time step
dt = 1e-9;
%Number of steps for timeout
nT = 3e3;
%Number of timesteps per ionization check
ionization_factor = 100;
%Ionization
file_inz = 'ADAS/scd93_c.dat';
%Recombination
file_rcmb = 'ADAS/acd93_c.dat';
