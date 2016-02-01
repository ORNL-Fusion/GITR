clear variables

constants
gimpInput

plotTracks = 1;

end_pos = dlmread('end_positions.txt',' ');  
run_param = dlmread('run_param.txt', ' ');
params = num2cell(run_param);
 [nS, dt] = deal(params{:});





Bfield3D.x = reshape(dlmread('Bfield_x_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.y = reshape(dlmread('Bfield_y_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.z = reshape(dlmread('Bfield_z_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.mag = reshape(dlmread('Bfield_mag_out.txt','\t'),nXv, nYv, nZv);
        
        temp_eV = reshape(dlmread('temp_eV_out.txt','\t'),nXv, nYv, nZv,nS);
        
        density_m3 = reshape(dlmread('density_m3_out.txt','\t'),nXv, nYv, nZv,nS);
        
        flowVelocity_ms.x = reshape(dlmread('flowVelocity_x_out.txt','\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.y = reshape(dlmread('flowVelocity_y_out.txt','\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.z = reshape(dlmread('flowVelocity_z_out.txt','\t'),nXv, nYv, nZv,nS);
        
        Efield3D.x = reshape(dlmread('Efield_x_out.txt','\t'),nXv, nYv, nZv);
        Efield3D.y = reshape(dlmread('Efield_y_out.txt','\t'),nXv, nYv, nZv);
        Efield3D.z = reshape(dlmread('Efield_z_out.txt','\t'),nXv, nYv, nZv);
        if plotTracks
        xHistory =  dlmread('xHistory_out.txt',' ');
        yHistory =  dlmread('yHistory_out.txt',' ');
        zHistory =  dlmread('zHistory_out.txt',' ');
        vxHistory =  dlmread('vxHistory_out.txt',' ');
        vyHistory =  dlmread('vyHistory_out.txt',' ');
        vzHistory =  dlmread('vzHistory_out.txt',' ');
        
        Z_History = dlmread('Z_History_out.txt', ' ');
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

impurityDensityTally = dlmread('impurityDensityTally_out.txt', '\t');


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