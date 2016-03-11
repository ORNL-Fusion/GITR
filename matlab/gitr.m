clear variables


% Load physical constants

constants

% Load run input file

gitrInput

disp(['Number of particles: ', num2str(nP)])
disp(['Number of time steps: ', num2str(nT)])

p = gcp('nocreate');

if isempty(p) == 1
% Create cluster objec
pc = parcluster('local');

% Make temporary folder for storing the cluster control files 
scratch_dir = 'scratch';
mkdir(scratch_dir);

% Explicitly set the JobStorageLocation to the temp directory that you
% have created in this script (above)
pc.JobStorageLocation = scratch_dir;

% Start the parallel pool with the number workers that matches your job
poolobj = parpool(pc);
end
% Load Processed ADAS data (cross sections, ionization rates, etc)

load('ADAS/IonizationData.mat');
load('ADAS/RecombinationData.mat');

% Create volume grid

yMinV = yMin;
yMaxV = yMax;

zMinV = zMin;
zMaxV = zMax;

xV_1D = linspace(xMinV,xMaxV,nXv);
yV_1D = linspace(yMinV,yMaxV,nYv);
zV_1D = linspace(zMinV,zMaxV,nZv);

dXv = xV_1D(2)-xV_1D(1);
dYv = yV_1D(2)-yV_1D(1);
dZv = zV_1D(2)-zV_1D(1);

xyz.x = xV_1D; % Create volume coordinate strucutre
xyz.y = yV_1D;
xyz.z = zV_1D;

% Create surface grid

surf_y1D = linspace(yMin,yMax,nY);
surf_z1D = linspace(zMin,zMax,nZ);
surf_x2D = zeros(nY,nZ);
surf_hist = zeros(nY,nZ);


for j=1:nY
    surf_x2D(j,:) =  (surf_z1D - surface_zIntercept) / surface_dz_dx;
end

% Plot initial (zeroed) surface histogram

if plotInitialSurface
    figure(1)
    h1 =  surf(surf_z1D,surf_y1D,surf_x2D,surf_hist);
    xlabel('z axis')
    ylabel('y axis')
    zlabel('x axis')
    title('Surface')
    axis equal
    drawnow
    xlim([zMinV zMaxV])
    ylim([yMinV yMaxV])
    zlim([xMinV xMaxV])
    az = 0;
    el = 0;
    view(az,el);
    set(gca,'Zdir','reverse')
    az = 45;
    el = 45;
    view(az,el);
end

% Setup B field

Bfield3D.x = zeros(nXv,nYv,nZv);
Bfield3D.y = zeros(nXv,nYv,nZv);
Bfield3D.z = zeros(nXv,nYv,nZv);
Bfield3D.mag = zeros(nXv,nYv,nZv);

BfieldInitialization

% Setup background plasma profiles - temperature, density, flow velocity

[n_, nS] = size(background_amu);

density_m3 = zeros(nXv,nYv,nZv,nS);
temp_eV = zeros(nXv,nYv,nZv,nS);
flowVelocity_ms.x = zeros(nXv,nYv,nZv,nS);
flowVelocity_ms.y = zeros(nXv,nYv,nZv,nS);
flowVelocity_ms.z = zeros(nXv,nYv,nZv,nS);

backgroundPlasmaInitialization


% Setup E field
debyeLength = sqrt(EPS0*maxTemp_eV/(maxDensity*background_Z(2)^2*Q));%only one q in order to convert to J

Efield3D.x = zeros(nXv,nYv,nZv);
Efield3D.y = zeros(nXv,nYv,nZv);
Efield3D.z = zeros(nXv,nYv,nZv);

EfieldInitialization


% Setup perpedicular diffusion coefficient
perDiffusionCoeff = zeros(nXv,nYv,nZv);

perDiffusionCoeffInitialization


% Plot slices through the profiles
if plot1DProfileSlices
    figure(2)
    subplot(2,2,1)
    semilogy(xV_1D,density_m3(:,1,1,1))
    title('Electron density')
    subplot(2,2,2)
    semilogy(xV_1D,temp_eV(:,1,1,1))
    title('Electron temp [eV]')
    subplot(2,2,3)
    plot(zV_1D,flowVelocity_ms.z(:,1,1,2))
    title('Ion flow velocity in z')
    subplot(2,2,4)
    plot(xV_1D,Efield3D.x(:,1,1))
    title('Ex [V/m]')
end

disp('Initialization of fields complete... Beginning particle initialization')
% Populate the impurity particle list
particles(nP) = particle;

for p=1:nP
    particles(p).Z = impurity_Z;
    particles(p).amu = impurity_amu;
    
    particles(p).x = x_start;
    particles(p).y = y_start;
    particles(p).z = z_start;
    
    particles(p).vx = sign(energy_eV_x_start) * sqrt(2*abs(energy_eV_x_start*Q)/(particles(p).amu*MI));
    particles(p).vy = sign(energy_eV_y_start) * sqrt(2*abs(energy_eV_y_start*Q)/(particles(p).amu*MI));
    particles(p).vz = sign(energy_eV_z_start) * sqrt(2*abs(energy_eV_z_start*Q)/(particles(p).amu*MI));
    
    particles(p).hitWall = 0;
    particles(p).leftVolume = 0;
    
    [s1,s2,s3,s4,s5,s6] = RandStream.create('mrg32k3a','NumStreams',6,'Seed','shuffle'); %Include ,'Seed','shuffle' to get different values each time
    
    particles(p).streams.ionization = s1;
    particles(p).streams.recombination = s2;
    particles(p).streams.perDiffusion = s3;
    particles(p).streams.parVelocityDiffusion = s4;
    particles(p).streams.per1VelocityDiffusion = s5;
    particles(p).streams.per2VelocityDiffusion = s6;
    
    particles(p).UpdatePrevious();
    
    particles(p).PerpDistanceToSurface(surface_dz_dx,surface_zIntercept);
    
end

% Calculate time step (dt)

% max_B = max( Bfield3D.mag(:) );
% min_m = impurity_amu * MI;
% max_wc = impurity_Z * Q * max_B / min_m;
% dt = 2 * pi / max_wc / nPtsPerGyroOrbit;
    dt = 1e-6/nPtsPerGyroOrbit;


% Setup arrays to store history

pre_history
nChargeStates = length(densityChargeBins);
nDensityBins = nChargeStates+1; %The chosen charge states and total
impurityDensityTally = zeros(nXv,nYv,nZv,nDensityBins);


% Main loop

IonizationTimeStep = ionization_nDtPerApply*dt;
CollisionTimeStep = collision_nDtPerApply*dt;

disp('Particle initialization complete... Starting main loop')

tic

% Progress bar
barLength = nP/progressInterval;
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,barLength) '\n\n']);

parfor p=1:nP
    
    tmp = zeros(nXv,nYv,nZv,nDensityBins);
    
    for tt = 1:nT
        
        particles(p).borisMove(xyz,Efield3D,Bfield3D,dt,...
            interpolators,positionStepTolerance,velocityChangeTolerance, ...
            debyeLength, -3*maxTemp_eV,surface_dz_dx,surface_zIntercept,background_Z,background_amu,maxTemp_eV);
        
        %                 [T Y] =  particles(p).move(dt,dt,Efield3D,Bfield3D,xyz,...
        %                     interpolators,xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV,...
        %                     surface_zIntercept,surface_dz_dx, ...
        %                     debyeLength, -3*maxTemp_eV,background_Z,background_amu,maxTemp_eV);
        
        particles(p).CrossFieldDiffusion(xyz,Bfield3D,perDiffusionCoeff,...
            interpolators,dt,positionStepTolerance);
        
        if particles(p).hitWall == 0 && particles(p).leftVolume ==0
            particles(p).UpdatePrevious();
            % particles(p).perpDistanceToSurface;
            
            [Mx, xIndex] = min(abs(xyz.x-particles(p).x));
            [My, yIndex] = min(abs(xyz.y-particles(p).y));
            [Mz, zIndex] = min(abs(xyz.z-particles(p).z));
            
            
            % This is to account for matlab not storing
            % structure arrays in parfor
            comp = find(particles(p).Z == densityChargeBins);
            if comp
                tmp(xIndex,yIndex,zIndex,comp) = tmp(xIndex,yIndex,zIndex,comp)+  dt;
            end
            tmp(xIndex,yIndex,zIndex,end) = tmp(xIndex,yIndex,zIndex,end)+  dt;
            
        end
        
        particles(p).OutOfDomainCheck(xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV);
        
%         particles(p).HitWallCheck(surface_zIntercept,surface_dz_dx);
         if  particles(p).hitWall == 0 && particles(p).leftVolume ==0
            
            particles(p).ionization(dt,xyz,density_m3,temp_eV,...
                IonizationData,interpolators,ionizationProbabilityTolerance);
            
            particles(p).recombination(dt,xyz,density_m3,temp_eV,...
                RecombinationData,interpolators,ionizationProbabilityTolerance);
        end
        
        
        
        if mod(tt, collision_nDtPerApply) == 0
            
            diagnostics = particles(p).CoulombCollisions(xyz,Bfield3D,flowVelocity_ms,density_m3,temp_eV,...
                background_amu,background_Z,interpolators, ...
                CollisionTimeStep,velocityChangeTolerance, connectionLength,surface_dz_dx,surface_zIntercept);
            
        end
        if trackHistory
            history(tt,p).x = particles(p).xPrevious;
            history(tt,p).y = particles(p).yPrevious;
            history(tt,p).z = particles(p).zPrevious;
            history(tt,p).vx = particles(p).vxPrevious;
            history(tt,p).vy = particles(p).vyPrevious;
            history(tt,p).vz = particles(p).vzPrevious;
            history(tt,p).Z = particles(p).Z;
        end
        
    end
    particles(p).Energy;
    % This is to account for matlab not storing
    % structure arrays in parfor - propagate worker changes back to client
    particles(p) = particles(p);
    impurityDensityTally = impurityDensityTally + tmp;
    
    if mod(p,progressInterval) == 0
        fprintf('\b|\n');
    end
end
toc


status = mkdir('output');

if trackHistory
    save('output/gitrHistories.mat','history');
end

run_param.nS = nS;
run_param.dt = dt;
run_param.surf_x2D = surf_x2D;
run_param.surf_y1D = surf_y1D;
run_param.surf_z1D = surf_z1D;
run_param.xV_1D = xV_1D;
run_param.yV_1D = yV_1D;
run_param.zV_1D = zV_1D;

save('output/gitrRunParameters.mat','run_param');
save('output/gitrParticles.mat','particles');
save('output/gitrImpurityDensityTally.mat','impurityDensityTally');

print_profiles
% At the completion of the job tidy up
% delete(poolobj);
% %delete(gcp('nocreate'))
% 
% % Clean up the temporary files
% rmdir(scratch_dir,'s');
%quit
