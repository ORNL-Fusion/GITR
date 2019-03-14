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

vGITR = [750 900 1000 1200 1400];
nGITR = 1e-4*[18882 19377 14996 11035 10133];
hold on
scatter(vGITR,nGITR)