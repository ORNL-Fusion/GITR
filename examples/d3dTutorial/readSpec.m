specFile = 'spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridZ = ncread(specFile,'gridZ');

total = dens(:,:,end);
pcolor(gridR,gridZ,total')
% colorbar
% set(gca,'YDir','normal')
% %axis equal
%  axis([0 0.2 0 .06])
% title('Radial Slice of Singly Ionized Tungsten Density [m-3]')
% xlabel('Z [m]')
% ylabel('R [m]')
% set(gca,'FontSize',15)