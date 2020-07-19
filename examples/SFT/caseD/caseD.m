close all
clear all
file = 'output/positions.nc';

hit = ncread(file,'hitWall');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');

figure(1)
histogram(vz)
hold on
histogram(vx)
histogram(vy)

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
% semilogy(s,boltzmann,'lineWidth',2)
title('Case D Impurity Density Profile')
xlabel('s [m]')
ylabel('C^{4+} density [m^{-3}]')
% axis([ 0 9 1e13 1e21])
set(gca,'fontsize',16)
% legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}')
dirs = {'outputDiamond','outputSquare','outputTriangle'};
dirs = {'output'};

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
boltzmann = 1e10*exp(-gridZ)
gitr_dens = 8e2*sum(dens(:,:,5),1);
semilogy(gridZ,boltzmann,'lineWidth',2)
semilogy(gridZ,gitr_dens)
axis([ 0 9 1e13 1e21])
set(gca, 'YScale', 'log')
end
legend('Simple Fluid Theory','T_D=10eV n=10^{19} v_D=-112m/s','T_D=100eV n=10^{19} v_D=-3.5e5e4m/s','T_D=100eV n=10^{18} v_D=-3.5e5e5m/s')

[v i] = max(gitr_dens);
figure(123)
scatter(gitr_dens(i:end),boltzmann(i:end))
b1 = gitr_dens(i:end)'\boltzmann(i:end);
hold on
plot(gitr_dens(i:end),gitr_dens(i:end)*b1)

yCalc1 = gitr_dens(i:end)*b1;
y = boltzmann(i:end)';
Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2)
function tau_s = drag_tau(mi,Zi,mD,TD,nD)
nD = nD/1e18;
lnGam = 15;
tau_s = mi*TD*(TD/mD)^(1/2)/(6.8e4*(1+mD/mi)*nD*Zi^2*lnGam);
end