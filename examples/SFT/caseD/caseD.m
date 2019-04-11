close all
clear all

TD = 10;
Z=4;
m=12;
sinj = 0.15;
n0 = 1e19;

    nn = n0
tau1 = drag_tau(12,4,2,TD,nn)
vTh = sqrt(TD*1.602e-19/12/1.66e-27);

D_par = vTh^2*tau1
s = linspace(0,10);
vD = -112.2;
vdiff = D_par/sinj;
vpl = vdiff+vD;
phi_in = 1.73e23;
np =phi_in/vpl; 
FFf = m*1.66e-27*(vD)./TD./1.602e-19./tau1.*(s);
ns = np*exp(FFf);
s = linspace(0,6,100);
figure(1)
hold on
semilogy(s,ns)

set(gca, 'YScale', 'log')

s = linspace(0,9,100);
figure(4)
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
semilogy(gridZ,8e12*sum(dens(:,:,5),1))
end
legend('Simple Fluid Theory','T_D=10eV n=10^{19} v_D=-112m/s','T_D=100eV n=10^{19} v_D=-3.5e5e4m/s','T_D=100eV n=10^{18} v_D=-3.5e5e5m/s')

function tau_s = drag_tau(mi,Zi,mD,TD,nD)
nD = nD/1e18;
lnGam = 15;
tau_s = mi*TD*(TD/mD)^(1/2)/(6.8e4*(1+mD/mi)*nD*Zi^2*lnGam);
end