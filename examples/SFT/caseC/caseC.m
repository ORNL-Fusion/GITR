close all
clear all

TD = 50;
E=-25;
Z=4;
m=12;
sinj = 0.15;
n0 = [1e18,1e19,1e20];
for i=1:length(n0)
    nn = n0(i);
tau1 = drag_tau(12,4,2,TD,nn)
vTh = sqrt(TD*1.602e-19/12/1.66e-27);
0.5*12*1.66e-27*vTh*vTh/1.602e-19
D_par = vTh^2*tau1
s = linspace(0,10);
vE = tau1*Z*E/m*1.602e-19/1.66e-27
vdiff = D_par/sinj;
vpl = vdiff+vE;
phi_in = 1.73e23;
np =phi_in/vpl; 
FEf = Z*E/TD*(s-sinj);
ns = np*FEf;
s = linspace(0,6,100);
figure(1)
hold on
semilogy(s,ns)
end
set(gca, 'YScale', 'log')
boltzmann = 3e19*exp(Z*E*s./TD)
figure(2)
semilogy(s,boltzmann,'lineWidth',2)
title('Case C Impurity Density Profile')
xlabel('s [m]')
ylabel('C^{4+} density [m^{-3}]')
% axis([0 150 0 1e20])
set(gca,'fontsize',16)
% legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}')
dirs = {'output1e18','output1e19','output1e20'};

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
semilogy(gridZ,2e13*sum(dens(:,:,5),1))
end
legend('Simple Fluid Theory','n_e = 10^{18} m^{-3}','n_e = 10^{19} m^{-3}','n_e = 10^{20} m^{-3}')
axis([0 6 1e13 1e20])
function tau_s = drag_tau(mi,Zi,mD,TD,nD)
nD = nD/1e18;
lnGam = 15;
tau_s = mi*TD*(TD/mD)^(1/2)/(6.8e4*(1+mD/mi)*nD*Zi^2*lnGam);
end