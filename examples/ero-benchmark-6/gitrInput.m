%input file for GITR0.0

% Impurity particles 

nP = 1e4;
sourceStrength = 1e19;

x_start = 0.00;
y_start = 0.00;
z_start = 0.00;

energy_eV_x_start = -10.0;
energy_eV_y_start = 0.0;
energy_eV_z_start = 0.0;

impurity_amu = 184.0;
impurity_Z = 0.0;

densityChargeBins = [0 1];

% Volume definition

xMinV =-0.002;
xMaxV =+0.002;

% Surface definition

yMin = -0.003;
yMax = 0.003;

zMin = -0.003;
zMax = 0.003;

% Volume grid
nXv = 100;
nYv = 150;
nZv = 150;

% Surface grid

nY = 150;
nZ = 150;

% Surface parameterization z = dz/dx * x + b

surface_dz_dx = 1.73205;
surface_zIntercept = 0;

connectionLength = 50;

% Background species info

background_Z = [-1 1];
background_amu = [ME/MI 2];

%Ionization
file_inz = 'ADAS/scd50_w.dat';
%Recombination
file_rcmb = 'ADAS/acd50_w.dat';
%Emission
file_emission = {'ADAS/w0_400875.m','ADAS/w1_434811.m'};

% Particle time stepping control

nPtsPerGyroOrbit = 50;
ionization_nDtPerApply = 1;
collision_nDtPerApply = 5;
nT = 5e3;

% Interpolator Dimensionality Selection

EfieldInterpolator_number = 2;
BfieldInterpolator_number = 0;
FlowVelocityInterpolator_number = 0;

temperatureInterpolator_number = 0;
densityInterpolator_number = 0;
perDiffusionCoefficientInterpolator_number = 0;

% Constant E field value - only used when EfieldInterpolator_number = 0
Efield_in = [1e2 0 0];

% Constant B field value - only used when BfieldInterpolator_number = 0
Bx_in = 0.00;
By_in = 0.00;
Bz_in = -2.0;

%Fraction of thermal velocity - only used when FlowVelocityInterpolator_number = 0
background_flow = [0.0 0.0];

% Perp DiffusionCoeff - only used when Diffusion interpolator is = 0
perDiffusionCoeff_in = 0.0;

% Background profile values used Density, temperature interpolators are 0
% or 2
maxDensity = [1e19 1e19];
densitySOLDecayLength =1e4;
maxTemp_eV = [20 20];
tempSOLDecayLength = 1e4;

% Checks on Monte Carlo Probability and Step Size
ionizationProbabilityTolerance = 2.0;
velocityChangeTolerance = 1e-2; % Fraction of previous speed
positionStepTolerance = 1e-3;

% Output options
printProfiles = 0;
trackHistory = 0;

% Plots

plotInitialSurface = 0;
plot1DProfileSlices = 0;

progressInterval = ceil(0.02*nP);
