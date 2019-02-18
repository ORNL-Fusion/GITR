clear all
close all
specFile = 'output/spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridZ = ncread(specFile,'gridZ');

total = dens(:,:,end);
figure(1)
p1=pcolor(gridR,gridZ,log10(total'))
set(p1, 'EdgeColor', 'none')
colorbar
% set(gca,'YDir','normal')
axis equal
%  axis([0 0.2 0 .06])
title({'ITER Ne Density [m-3]','Other parameters'})
xlabel('Z [m]')
ylabel('R [m]')
set(gca,'FontSize',15)
M = csvread('iterGeom.csv');
r = M(:,1);
z = M(:,2);
hold on
plot([r;r(1)],[z;z(1)],'w','lineWidth',2)
axis([min(gridR) max(gridR) min(gridZ) max(gridZ)])