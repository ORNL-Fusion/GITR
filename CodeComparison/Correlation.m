axisLabelFont = 20;
tickFont = 12;
titleFont = 60;

dA = CellProp{6}.data;
Come = NC*dA(1,1);
R = corrcoef(rot90(rot90(Tallys)),Come)

figure(11)
gitrDeposition = rot90(rot90(Tallys));
set(gcf,'Position',[300 300 900 600])
subplot(2,1,1)
surf(x,y,zeros(iNX,iNY),log10(gitrDeposition),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -5 8])
caxis([-1 3.5])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


subplot(2,1,2)
surf(x,y,zeros(iNX,iNY),log10(Come),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -5 8])
caxis([-1 3.5])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


% subplot(2,2,2)
% surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_MQ)),'EdgeColor','none')
% colormap hot
% colorbar
% axis([-2 15 -5 8])
% caxis([0 6])
% xlabel('x axis [mm]','FontSize',axisLabelFont)
% ylabel('y axis [mm]','FontSize',axisLabelFont)
% zlabel('z axis [mm]','FontSize',axisLabelFont)
% title('GITR Mean Charge [#]','FontSize',titleFont)
% set(gca,'FontSize',tickFont)
% 
% subplot(2,2,4)
% surf(x,y,zeros(iNX,iNY),MQ,'EdgeColor','none')
% colormap hot
% colorbar
% axis([-2 15 -5 8])
% caxis([0 6])
% xlabel('x axis [mm]','FontSize',axisLabelFont)
% ylabel('y axis [mm]','FontSize',axisLabelFont)
% zlabel('z axis [mm]','FontSize',axisLabelFont)
% title('ERO Mean Charge [#]','FontSize',titleFont)
% set(gca,'FontSize',tickFont)

% subplot(2,3,3)
% surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_ME)),'EdgeColor','none')
% colormap hot
% colorbar
% axis([-2 15 -5 8])
% caxis([0 350])
% xlabel('x axis [mm]','FontSize',axisLabelFont)
% ylabel('y axis [mm]','FontSize',axisLabelFont)
% zlabel('z axis [mm]','FontSize',axisLabelFont)
% title('GITR Mean Impact Energy [eV]','FontSize',titleFont)
% set(gca,'FontSize',tickFont)
% 
% subplot(2,3,6)
% surf(x,y,zeros(iNX,iNY),ME,'EdgeColor','none')
% colormap hot
% colorbar
% axis([-2 15 -5 8])
% caxis([0 350])
% xlabel('x axis [mm]','FontSize',axisLabelFont)
% ylabel('y axis [mm]','FontSize',axisLabelFont)
% zlabel('z axis [mm]','FontSize',axisLabelFont)
% title('ERO Mean Impact Energy [eV]','FontSize',titleFont)
% set(gca,'FontSize',tickFont)

whitebg('white')
figure(12)
%histogram2(Come,gitrDeposition)

xx = reshape(gitrDeposition,[150*150,1]);
yy = reshape(Come, [150*150, 1]);

b1 = xx\yy;
scatter(log10(xx),log10(yy))

vals = 0:1:700;
bs = vals*b1;
hold on
plot(log10(vals),log10(bs))
xlabel('log10 GITR Deposition Values','FontSize',axisLabelFont)
ylabel('log10 ERO Deposition Values','FontSize',axisLabelFont)
title('Least Squares Fit for GITR-ERO Comparison of Deposition','FontSize',titleFont)
set(gca,'FontSize',tickFont)
hold off
yCalc1 = b1*xx;
Rsq1 = 1 - sum((yy - yCalc1).^2)/sum((yy - mean(yy)).^2)