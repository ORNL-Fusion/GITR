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