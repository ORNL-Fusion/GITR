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
Bx_in = +0.1;
% Bfieldy
By_in = +1.0;
% Bfieldz
Bz_in = -0.3;

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

x_start = -0.004;
y_start = -0.01;
z_start = 0.0;

energy_eV_x_start = 1.0;
energy_eV_y_start = 1.5;
energy_eV_z_start = -1.0;

impurity_amu = 12.0;
impurity_Z = 1.0;

%Slope of parameterized surface (X = m*Y +b)
surf_slope = -0.1;
%Intercept of parameterized surface
surf_incpt = -0.002;


%Ionization
file_inz = 'ADAS/scd93_c.dat';
%Recombination
file_rcmb = 'ADAS/acd93_c.dat';

% Particle time stepping control

nPtsPerGyroOrbit = 1e3;
ionization_nDtPerApply = 100;
%dt = 1e-9;
max_nT = 3000;

% Plots

plotInitialSurface = 1;
plot1DProfileSlices = 1;
