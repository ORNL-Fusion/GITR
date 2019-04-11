clear all
close all
caseAanalysis
vGITR = [750,900,1000,1200,1400]
close all
vB = linspace(700,1600);
np = 1.73e23./(243+vB);

figure(1)
plot(vB,np./1e20,'lineWidth',2)
title('Case A Peak Impurity Density')
xlabel('-v_B [m/s]')
ylabel('n_p [10^{20} m^{-3}]')
axis([700 1600 0 2])
set(gca,'fontsize',16)
hold on
scatter(vGITR,npGITR/nGITR_noScatt(1)*np(7)/1e20,'k')

% M = csvread('../SFT-2/caseA-2-Table 1.csv',1,0)
% vGITR = M(:,1);
% nGITR_noScatt = M(:,6);
% nGITR_Scatt = M(:,7);
% hold on
% scatter(vGITR,nGITR_noScatt/nGITR_noScatt(1)*np(7)/1e20,'k')
% s=scatter(vGITR,nGITR_Scatt/nGITR_noScatt(1)*np(7)/1e20,'d')
% s.MarkerEdgeColor = 'k';
% legend('Simple Fluid Theory','GITR - Drag and Heating', 'GITR - Drag, Heating, Pitch-Angle Scattering')
legend('Simple Fluid Theory','GITR')
pbaspect([1 1.25 1])