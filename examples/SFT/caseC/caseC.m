close all
clear all

TD = 50;
E=-25;
Z=4;
s = linspace(0,6,100);

boltzmann = 3e19*exp(Z*E*s./TD)
semilogy(s,boltzmann,'lineWidth',2)
title('Case C Impurity Density Profile')
xlabel('s [m]')
ylabel('C^{4+} density [m^{-3}]')
% axis([0 150 0 1e20])
set(gca,'fontsize',16)
% legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}')