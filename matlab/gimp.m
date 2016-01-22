clear variables


% Load physical constants

constants

% Load run input file

gimpInput

% Load ADAS data (cross sections, ionization rates, etc)

[IonizationTemp, IonizationDensity, IonizationRateCoeff, IonizationChargeState] = ADF11(file_inz);
[RecombinationTemp, RecombinationDensity, RecombinationRateCoeff, RecombinationChargeState] = ADF11(file_rcmb);

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

Bx = zeros(nXv,nYv,nZv);
By = zeros(nXv,nYv,nZv);
Bz = zeros(nXv,nYv,nZv);

Bx(:) = Bx_in;
By(:) = By_in;
Bz(:) = Bz_in;

Bfield3D.x = Bx; % Create background B strucutre
Bfield3D.y = By;
Bfield3D.z = Bz;

BMag = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Setup profiles

[n_, nS] = size(background_amu);

density_m3 = zeros(nXv,nYv,nZv,nS);
temp_eV = zeros(nXv,nYv,nZv,nS);

impurityDensity = zeros(nT,nP);
volumeGridSize = [nXv nYv nZv];

Ex = zeros(nXv,nYv,nZv);
Ey = zeros(nXv,nYv,nZv);
Ez = zeros(nXv,nYv,nZv);

Efield3D.x = Ex; % Create E field structure
Efield3D.y = Ey;
Efield3D.z = Ez;

perDiffusionCoeff = zeros(nXv,nYv,nZv);

V_1D = sheathPotential * exp( xV_1D / sheathWidth );

debyeLength = sqrt(EPS0*maxTemp_eV/(maxDensity*background_Z(2)^2*Q));%only one q in order to convert to J

for i=1:nXv
    for j=1:nYv
        for k=1:nZv
            
            density_m3(i,j,k,:) = maxDensity * exp( (xMinV-xV_1D(i)) / densitySOLDecayLength );
            temp_eV(i,j,k,:) = maxTemp_eV * exp( (xMinV-xV_1D(i)) / tempSOLDecayLength );
            
            perDiffusionCoeff(i,j,k) = perDiffusionCoeff_in;
            
            if i>1 && i<nXv
                Efield3D.x(i,j,k) = -(V_1D(i+1)-V_1D(i-1)) / (2*dXv);
            elseif i==1
                Efield3D.x(i,j,k) = -(-1*V_1D(i)+V_1D(i+1)) / dXv;
            elseif i==nXv
                Efield3D.x(i,j,k) = -(-V_1D(i-1)+V_1D(i)) / dXv;
            end
            
            Efield3D.y(i,j,k) = 0;
            Efield3D.z(i,j,k) = 0;
            
        end
    end
end


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
    plot(xV_1D,V_1D(:))
    title('Sheath potential [V]')
    subplot(2,2,4)
    plot(xV_1D,Efield3D.x(:,1,1))
    title('Ex [V/m]')
end

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
    
end

% Calculate time step (dt)

max_B = max( BMag(:) );
min_m = impurity_amu * MI;
max_wc = impurity_Z * Q * max_B / min_m;
dt = 2 * pi / max_wc / nPtsPerGyroOrbit;
if impurity_Z == 0
    dt = 1e-6/nPtsPerGyroOrbit;
end

% Setup arrays to store history

pre_history

% Main loop

IonizationTimeStep = ionization_nDtPerApply*dt;

tic
parfor p=1:nP
    
    for tt = 1:nT
        
        particles(p).PerpDistanceToSurface(surface_dz_dx,surface_zIntercept);
      

        if mod(tt, ionization_nDtPerApply) == 0
            
            particles(p).ionization(IonizationTimeStep,xyz,density_m3,temp_eV,...
               IonizationRateCoeff,IonizationTemp, IonizationDensity,...
               IonizationChargeState,selectedScalarInterpolator,ionizationProbabilityTolerance);
            
            particles(p).recombination(IonizationTimeStep,xyz,density_m3,temp_eV,...
               RecombinationRateCoeff,RecombinationTemp,RecombinationDensity,...
               RecombinationChargeState,selectedScalarInterpolator,ionizationProbabilityTolerance);
            
        end

%         particles(p).CrossFieldDiffusion(Bfield3D,xyz,perDiffusionCoeff,dt,...
%             selectedScalarInterpolator,selectedVectorInterpolator,positionStepTolerance);


        diagnostics = particles(p).CoulombCollisions(xyz,Bfield3D,density_m3,temp_eV,...
            background_amu,background_Z,dt,...
            selectedVectorInterpolator,selectedScalarInterpolator, ...
            velocityChangeTolerance, connectionLength,maxTemp_eV,surface_dz_dx,surface_zIntercept);

        particles(p).borisMove(xyz,Efield3D,Bfield3D,dt,...
            selectedVectorInterpolator,positionStepTolerance,velocityChangeTolerance, ...
            debyeLength, -3*maxTemp_eV,surface_dz_dx,background_Z,background_amu,maxTemp_eV);

%         [T Y] =  particles(p).move(dt,dt,Efield3D,Bfield3D,xyz,...
%             selectedVectorInterpolator,xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV,...
%             surface_zIntercept,surface_dz_dx)
        

        


        particles(p).OutOfDomainCheck(xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV);
        
        particles(p).HitWallCheck(surface_zIntercept,surface_dz_dx);
        

        if particles(p).hitWall == 0 && particles(p).leftVolume ==0
            particles(p).UpdatePrevious();
            
            [Mx, xIndex] = min(abs(xyz.x-particles(p).x));
            [My, yIndex] = min(abs(xyz.y-particles(p).y));
            [Mz, zIndex] = min(abs(xyz.z-particles(p).z));
             
             impurityDensity(tt,p) = sub2ind(volumeGridSize,xIndex,yIndex,zIndex);
        end
        
        
        xHistory(tt,p) = particles(p).xPrevious;
        yHistory(tt,p) = particles(p).yPrevious;
        zHistory(tt,p) = particles(p).zPrevious;
        vxHistory(tt,p) = particles(p).vxPrevious;
        vyHistory(tt,p) = particles(p).vyPrevious;
        vzHistory(tt,p) = particles(p).vzPrevious;
          
    end
    end_pos(p,:) = [particles(p).xPrevious particles(p).yPrevious particles(p).zPrevious];
end
toc

history_plot
surface_scatter

fileID = fopen('end_positions.txt','w');

fprintf(fileID,'%f %f %f\n',end_pos);
fclose(fileID);