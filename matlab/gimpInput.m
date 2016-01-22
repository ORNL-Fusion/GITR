%input file for GITR0.0

% Volume definition

xMinV =-0.05;
xMaxV =+0.05;

% Surface definition

yMin = -0.05;
yMax = 0.05;

zMin = -0.05;
zMax = 0.05;

% Surface grid

nY = 100;
nZ = 100;

% Volume grid
nXv = 50;
nYv = 50;
nZv = 50;

% Sheath potential
sheathPotential = -60.0;
% Potential decay length
sheathWidth = 0.0001;

Bx_in = 0.00;
By_in = 0.00;
Bz_in = -2.0;

% Background species info

background_Z = [-1 1];
background_amu = [ME/MI 2];

% Density
maxDensity = 1e19;
% Density Decay length
densitySOLDecayLength =1e4;
%Temperature (in eV)
maxTemp_eV = 20;
%Temperature decay length
tempSOLDecayLength = 1e4;
%Dperp
perDiffusionCoeff_in = 0.04;


% Impurity particles 

nP = 10;

x_start = 0.00;
y_start = 0.00;
z_start = 0.00;

energy_eV_x_start = -10.0;
energy_eV_y_start = 0.0;
energy_eV_z_start = 0.0;


impurity_amu = 184.0;
impurity_Z = 0.0;

% Surface parameterization z = dz/dx * x + b

surface_dz_dx = 5;
surface_zIntercept = 0;

%Ionization
file_inz = 'ADAS/scd50_w.dat';
%Recombination
file_rcmb = 'ADAS/acd50_w.dat';

% Particle time stepping control

nPtsPerGyroOrbit = 100;
ionization_nDtPerApply = 50;
nT = 2000;

% Plots

plotInitialSurface = 1;
plot1DProfileSlices = 1;

% Interpolator Dimensionality Selection

selectedVectorInterpolator = @gimpInterpVector1D;
selectedScalarInterpolator = @gimpInterpScalar1D;

%selectedInterpolator = @gimpInterVectorp3D;
%selectedInterpolator = @gimpInterpScalar3D;

% Checks on Monte Carlo Probability and Step Size
ionizationProbabilityTolerance = 0.9;
velocityChangeTolerance = 1; % Fraction of previous speed
positionStepTolerance = 1e-3;

connectionLength = 50;

