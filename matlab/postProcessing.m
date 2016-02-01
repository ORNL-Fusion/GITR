clear variables
close all
constants

HistoryFileName = 'output/gitrHistories.mat';
InputFileName = 'gitrInput.m'
RunParametersFileName = 'output/gitrRunParameters.mat';
ParticlesFileName = 'output/gitrParticles.mat';

if exist(InputFileName)
   run(InputFileName);
else
    error('Could not find GITR input file ... WTF?');
end

plotScatter = 0;
plotTracks = 0;
plotDens = 0;
plotProfiles = 0;
plot3D = 0;

load(RunParametersFileName);
load(ParticlesFileName);

if plotProfiles
Bfield3D.x = reshape(dlmread(full_path{3},'\t'),nXv, nYv, nZv);
        Bfield3D.y = reshape(dlmread(full_path{4},'\t'),nXv, nYv, nZv);
        Bfield3D.z = reshape(dlmread(full_path{5},'\t'),nXv, nYv, nZv);
        Bfield3D.mag = reshape(dlmread(full_path{6},'\t'),nXv, nYv, nZv);
        
        temp_eV = reshape(dlmread(full_path{7},'\t'),nXv, nYv, nZv,nS);
        
        density_m3 = reshape(dlmread(full_path{8},'\t'),nXv, nYv, nZv,nS);
        
        flowVelocity_ms.x = reshape(dlmread(full_path{9},'\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.y = reshape(dlmread(full_path{10},'\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.z = reshape(dlmread(full_path{11},'\t'),nXv, nYv, nZv,nS);
        
        Efield3D.x = reshape(dlmread(full_path{12},'\t'),nXv, nYv, nZv);
        Efield3D.y = reshape(dlmread(full_path{13},'\t'),nXv, nYv, nZv);
        Efield3D.z = reshape(dlmread(full_path{14},'\t'),nXv, nYv, nZv);
end
        if plotTracks

        end
        yMinV = yMin;
yMaxV = yMax;

zMinV = zMin;
zMaxV = zMax;

xV_1D = linspace(xMinV,xMaxV,nXv);
yV_1D = linspace(yMinV,yMaxV,nYv);
zV_1D = linspace(zMinV,zMaxV,nZv);
surf_y1D = linspace(yMin,yMax,nY);
surf_z1D = linspace(zMin,zMax,nZ);
surf_x2D = zeros(nY,nZ);
surf_hist = zeros(nY,nZ);

if plotDens
impurityDensityTally = dlmread(full_path{22}, '\t');
end

for j=1:nY
    surf_x2D(j,:) =  (surf_z1D - surface_zIntercept) / surface_dz_dx;
end

% Plot initial (zeroed) surface histogram


    figure(1)
    h1 =  surf(surf_z1D,surf_y1D,surf_x2D,surf_hist);
    xlabel('z axis')
    ylabel('y axis')
    zlabel('x axis')
    title('Surface')
    axis equal
    drawnow
    xlim([zMin zMax])
    ylim([yMin yMax])
    zlim([xMinV xMaxV])
    az = 0;
    el = 0;
    view(az,el);
    set(gca,'Zdir','reverse')
    az = 45;
    el = 45;
    view(az,el);


if exist(HistoryFileName)
    load(HistoryFileName);
    history_plot
end

% Surface particle impact histogram

histogram2(particles.y,particles.z,'DisplayStyle','tile','ShowEmptyBins','on')

surface_scatter

density_calc