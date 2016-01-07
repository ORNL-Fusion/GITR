%input file for GITR0.0

% Volume definition

xMinV =-0.005;
xMaxV =+0.005;

% Surface definition

yMin = -0.005;
yMax = 0.005;

zMin = -0.005;
zMax = 0.005;

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
Bx_in = +0.4;
% Bfieldy
By_in = +0.1;
% Bfieldz
Bz_in = -1.0;

% Background species info

background_Z = [-1 1];
background_amu = [ME/MI 2];

% Density
maxDensity = 1e19;
% Density Decay length
densitySOLDecayLength =0.1;
%Temperature (in eV)
maxTemp_eV = 20;
%Temperature decay length
tempSOLDecayLength = .1;
%Dperp
perDiffusionCoeff_in = 0.04;


% Impurity particles 

nP = 2;

x_start = -0.002;
y_start = 0.00;
z_start = 0.003;

energy_eV_x_start = 1.0;
energy_eV_y_start = 1.5;
energy_eV_z_start = -1.0;

impurity_amu = 12.0;
impurity_Z = 1.0;

% Surface parameterization z = dz/dx * x + b

surface_dz_dx = 5;
surface_zIntercept = 0;

%Ionization
file_inz = 'ADAS/scd93_c.dat';
%Recombination
file_rcmb = 'ADAS/acd93_c.dat';

% Particle time stepping control

nPtsPerGyroOrbit = 1e3;
ionization_nDtPerApply = 100;
nT = 3000;

% Plots

plotInitialSurface = 1;
plot1DProfileSlices = 1;
