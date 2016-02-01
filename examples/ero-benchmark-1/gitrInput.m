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
sheathWidth = 0.00001;
% Constant E field value - only used when EfieldInterpolator_number = 0
Efield_in = [1e2 0 0];

Bx_in = 0.00;
By_in = 0.00;
Bz_in = -2.0;

% Background species info

background_Z = [-1 1];
background_amu = [ME/MI 2];
background_flow = [0.99 0.99];%Fraction of thermal velocity - only used when FlowVelocityInterpolator_number = 0

% Density
maxDensity = [1e19 1e19];
% Density Decay length
densitySOLDecayLength =1e4;
%Temperature (in eV)
maxTemp_eV = [20 20];
%Temperature decay length
tempSOLDecayLength = 1e4;
%Dperp
perDiffusionCoeff_in = 0.000004;


% Impurity particles 

nP = 100;

x_start = 0.00;
y_start = 0.00;
z_start = 0.00;

energy_eV_x_start = -10.0;
energy_eV_y_start = 0.0;
energy_eV_z_start = 0.0;


impurity_amu = 184.0;
impurity_Z = 0.0;

% Surface parameterization z = dz/dx * x + b

surface_dz_dx = 1.73205;
surface_zIntercept = 0;

%Ionization
file_inz = 'ADAS/scd50_w.dat';
%Recombination
file_rcmb = 'ADAS/acd50_w.dat';

% Particle time stepping control

nPtsPerGyroOrbit = 1e3;
ionization_nDtPerApply = 1;
nT = 1e4;
sheath_timestep_factor = 1e4;

% Plots

plotInitialSurface = 0;
plot1DProfileSlices = 0;

trackHistory = 1;

% Interpolator Dimensionality Selection

EfieldInterpolator_number = 2;
BfieldInterpolator_number = 0;
FlowVelocityInterpolator_number = 0;

temperatureInterpolator_number = 0;
densityInterpolator_number = 0;
perDiffusionCoefficientInterpolator_number = 0;


% Checks on Monte Carlo Probability and Step Size
ionizationProbabilityTolerance = 0.5;
velocityChangeTolerance = 1; % Fraction of previous speed
positionStepTolerance = 1e-3;

connectionLength = 50;

% Output options
printProfiles = 0;
printHistory = 1;

    




    
