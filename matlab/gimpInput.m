%input file for GITR0.0

% Volume definition

xMinV =-0.005;
xMaxV =+0.005;

% Surface definition

yMin = -0.005;
yMax = 0.005;

zMin = -0.005;
zMax = 0.005;

% Surface grid

nY = 30;
nZ = 40;

% Volume grid
nXv = 50;
nYv = 80;
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

nP = 96;

x_start = xMinV*0.95;
y_start = 0.00;
z_start = 0.00;

energy_eV_x_start = 1.0;
energy_eV_y_start = 0;
energy_eV_z_start = 0;

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

nPtsPerGyroOrbit = 1e2;
ionization_nDtPerApply = 100;
nT = 200;

% Plots

plotInitialSurface = 1;
plot1DProfileSlices = 1;

% Interpolator Dimensionality Selection

selectedVectorInterpolator = @gimpInterpVector1D;
selectedScalarInterpolator = @gimpInterpScalar1D;

%selectedInterpolator = @gimpInterp3D;
