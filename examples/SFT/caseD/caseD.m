close all
clear all

s = linspace(0,9,100);

boltzmann = 3e19*exp(-s)
semilogy(s,boltzmann,'lineWidth',2)
title('Case D Impurity Density Profile')
xlabel('s [m]')
ylabel('C^{4+} density [m^{-3}]')
axis([ 0 9 1e13 1e21])
set(gca,'fontsize',16)
% legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}')
dirs = {'outputDiamond','outputSquare','outputTriangle'};

for i=1:length(dirs)
specFile = [dirs{i},'/spec.nc'];
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
% gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');
% figure(2)
% h=pcolor(gridR,gridZ,dens(:,:,5)')
% h.EdgeColor = 'none';
% figure(21)
hold on
semilogy(gridZ,4e12*sum(dens(:,:,5),1))
end
legend('Simple Fluid Theory','Diamond','Square','Triangle')