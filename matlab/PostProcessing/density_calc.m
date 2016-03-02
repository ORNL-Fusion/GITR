nChargeStates = length(densityChargeBins);
n = 4;

    dy = yV_1D(2) - yV_1D(1);
    dx = xV_1D(2) - xV_1D(1);
    dz = zV_1D(2) - zV_1D(1);
density = impurityDensityTally(:,:,:,1)*sourceStrength/(dx*dy*dz)/nP;
run(file_emission{1})
radiance = zeros(nXv,nYv,nZv);

[a,b,c] =  ind2sub(size(density),find(density ~= 0));
len1 = length(a);

for i=1:len1
%             density_m3(i,j,k,1);
%             temp_eV(i,j,k,1);
            Coeff = interpn(Ne,Te,vSig_w0,density_m3(a(i),b(i),c(i),1)/1e6,temp_eV(a(i),b(i),c(i),1),'linear',0);
            radiance(a(i),b(i),c(i)) = Coeff*density(a(i),b(i),c(i))/1e6*density_m3(a(i),b(i),c(i),1);
                
        end

%for i=1:nChargeStates

sightLine = sum(radiance,2)*dy;
sightLine = reshape(sightLine,[nXv,nZv]);

figure(n+1)
surf( run_param.zV_1D,run_param.xV_1D,sightLine)
    xlabel('z axis')
    ylabel('x axis')
    zlabel('y axis')
    number = num2str(1);
    title0 = strcat('Tungsten  ', number, ' Line intensity (4009 Angstrom) [Photons/sr/s/m^2]');
    title(title0)
        axis([-0.008 0.008 -0.013 0.002])
        colormap jet
% h2 = histogram2(density(:,1),density(:,2),surf_y1D,surf_z1D,...
%     'DisplayStyle','tile','ShowEmptyBins','on')
% colorbar

% xslice = 0;
% yslice = 0.0;
% zslice = 0;
% figure(4)
% slice(xV_1D,yV_1D,zV_1D,density,xslice,yslice,zslice)
% colorbar
% x1 =  (zV_1D - surface_zIntercept) / surface_dz_dx;
% y1 = yV_1D*0.0 + yslice;
 hold on
 plot3(run_param.surf_z1D,run_param.surf_x2D(1,:),0*run_param.surf_z1D)
 hold off
     az = 180;
    el = 90;
    view(az,el);
    colorbar

%end
%  
%      xlabel('x axis')
%     ylabel('y axis')
%     zlabel('z axis')
%     title('Surface')
%     axis equal




density = impurityDensityTally(:,:,:,2)*sourceStrength/(dx*dy*dz)/nP;
run(file_emission{2})
radiance = zeros(nXv,nYv,nZv);

[a,b,c] =  ind2sub(size(density),find(density ~= 0));
len1 = length(a);

for i=1:len1
%             density_m3(i,j,k,1);
%             temp_eV(i,j,k,1);
            Coeff = interpn(Ne,Te,vSig_w1,density_m3(a(i),b(i),c(i),1)/1e6,temp_eV(a(i),b(i),c(i),1),'linear',0);
            radiance(a(i),b(i),c(i)) = Coeff*density(a(i),b(i),c(i))/1e6*density_m3(a(i),b(i),c(i),1);
                
        end



sightLine = sum(radiance,2)*dy;
sightLine = reshape(sightLine,[nXv,nZv]);
figure(n+2)
surf(run_param.zV_1D, run_param.xV_1D,sightLine)
    xlabel('z axis')
    ylabel('x axis')
    zlabel('y axis')
    number = num2str(2);
    title0 = strcat('Tungsten  ', number, ' Line intensity (4348 Angstrom) [Photons/sr/s/m^2]');
    title(title0)
            axis([-0.01 0.010 -0.015 0.005])
            colormap jet
     hold on
 plot3(run_param.surf_z1D,run_param.surf_x2D(1,:),0*run_param.surf_z1D)
 hold off
     az = 180;
    el = 90;
    view(az,el);

    colorbar