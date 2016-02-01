clear variables
close all
constants
% filePath = '~/GitHub/gimp/matlab/data/01272016/'
% inputFile = strcat(filePath,'gimpInput.m')
%  run inputFile


[FileName,PathName,FilterIndex] = uigetfile;
inpPath = strcat(PathName,FileName);
run(inpPath)
plotScatter = 1;
plotTracks = 1;
plotDens = 0;
plotProfiles = 1;
filenames = {'end_positions.txt','run_param.txt','Bfield_x_out.txt','Bfield_y_out.txt','Bfield_z_out.txt','Bfield_mag_out.txt', ...
   'temp_eV_out.txt','density_m3_out.txt','flowVelocity_x_out.txt','flowVelocity_y_out.txt','flowVelocity_z_out.txt', ...
   'Efield_x_out.txt','Efield_y_out.txt','Efield_z_out.txt','xHistory_out.txt','yHistory_out.txt','zHistory_out.txt', ...
   'vxHistory_out.txt','vyHistory_out.txt','vzHistory_out.txt','Z_History_out.txt','impurityDensityTally_out.txt'};
for i=1:numel(filenames)
full_path{i} = fullfile ( PathName, filenames{i} );
end
end_pos = dlmread(full_path{1},' ');  
run_param = dlmread(full_path{2}, ' ');
params = num2cell(run_param);
 [nS, dt] = deal(params{:});




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
        xHistory =  dlmread(full_path{15},' ');
        yHistory =  dlmread(full_path{16},' ');
        zHistory =  dlmread(full_path{17},' ');
        vxHistory =  dlmread(full_path{18},' ');
        vyHistory =  dlmread(full_path{19},' ');
        vzHistory =  dlmread(full_path{20},' ');
        
        Z_History = dlmread(full_path{21}, ' ');
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


history_plot
surface_scatter

density_calc