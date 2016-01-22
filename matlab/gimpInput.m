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

nPtsPerGyroOrbit = 70;
ionization_nDtPerApply = 50;
nT = 1500;

% Plots

plotInitialSurface = 1;
plot1DProfileSlices = 1;

% Interpolator Dimensionality Selection

EfieldInterpolator = 2;
BfieldInterpolator = 0;
FlowVelocityInterpolator = 2;

temperatureInterpolator = 0;
densityInterpolator = 0;


% Checks on Monte Carlo Probability and Step Size
ionizationProbabilityTolerance = 0.9;
velocityChangeTolerance = 1; % Fraction of previous speed
positionStepTolerance = 1e-3;

connectionLength = 50;

if EfieldInterpolator ==0
    EfieldInterpolator = @gitrEfieldConstant;
else if EfieldInterpolator == 1
        EfieldInterpolator = @gitrInterpVector1D;
    else if EfieldInterpolator == 2
            EfieldInterpolator = @gitrEfieldAnalytic;
        else if EfieldInterpolator == 3
                EfieldInterpolator = @gimpInterVectorp3D;
            end
        end
    end
end

if BfieldInterpolator ==0
    BfieldInterpolator = @gitrBfieldConstant;
else if EfieldInterpolator == 1
        BfieldInterpolator = @gitrInterpVector1D;
    else if BfieldInterpolator == 2
            BfieldInterpolator = @gitrBfieldAnalytic;
        else if BfieldInterpolator == 3
                BfieldInterpolator = @gimpInterVectorp3D;
            end
        end
    end
end
    
if FlowVelocityInterpolator ==0
    FlowVelocityInterpolator = @gitrFlowVelocityConstant;
else if FlowVelocityInterpolator == 1
        FlowVelocityInterpolator = @gitrInterpVector1D;
    else if FlowVelocityInterpolator == 2
            FlowVelocityInterpolator = @gitrFlowVelocityAnalytic;
        else if FlowVelocityInterpolator == 3
                FlowVelocityInterpolator = @gimpInterVectorp3D;
            end
        end
    end
end
    
if temperatureInterpolator ==0
    temperatureInterpolator = @gitrTemperatureConstant;
else if EfieldInterpolator == 1
        temperatureInterpolator = @gitrInterpScalar1D;
    else if temperatureInterpolator == 2
            temperatureInterpolator = @gitrTemperatureAnalytic;
        else if temperatureInterpolator == 3
                temperatureInterpolator = @gimpInterpScalar3D;
            end
        end
    end
end

if densityInterpolator ==0
    densityInterpolator = @gitrDensityConstant;
else if EfieldInterpolator == 1
        densityInterpolator = @gitrInterpScalar1D;
    else if densityInterpolator == 2
            densityInterpolator = @gitrDensityAnalytic;
        else if densityInterpolator == 3
                densityInterpolator = @gimpInterpScalar3D;
            end
        end
    end
end

    
interpolators = {EfieldInterpolator; BfieldInterpolator;...
    FlowVelocityInterpolator; temperatureInterpolator; ...
    densityInterpolator};