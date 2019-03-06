clear all
close all

vB = linspace(700,1600);
np = 1.73e23./(243+vB);

figure(1)
plot(vB,np./1e20,'lineWidth',2)
title('Case A Peak Impurity Density')
xlabel('-v_B [m/s]')
ylabel('n_p [10^{20} m^{-3}]')
axis([700 1600 0.2 2])
set(gca,'fontsize',16)