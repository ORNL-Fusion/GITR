axisLabelFont = 20;
tickFont = 12;
titleFont = 60;

dA = CellProp{6}.data;
Come = NC*dA(1,1);
R = corrcoef(rot90(rot90(Tallys)),Come)

figure(11)
gitrDeposition = rot90(rot90(Tallys));
set(gcf,'Position',[300 300 500 600])
subplot(2,3,1)
surf(x,y,zeros(iNX,iNY),log10(gitrDeposition),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([-1 3])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


subplot(2,3,4)
surf(x,y,zeros(iNX,iNY),log10(Come),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([-1 3])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


subplot(2,3,2)
surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_MQ)),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([0 6])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Mean Charge [#]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,5)
surf(x,y,zeros(iNX,iNY),MQ,'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([0 6])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Mean Charge [#]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,3)
surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_ME)),'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([0 350])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Mean Impact Energy [eV]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,6)
surf(x,y,zeros(iNX,iNY),ME,'EdgeColor','none')
colormap hot
colorbar
axis([-2 15 -2 8])
caxis([0 350])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Mean Impact Energy [eV]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

whitebg('white')
figure(12)
%histogram2(Come,gitrDeposition)
sizes  = length(x)*length(y)
xx = reshape(gitrDeposition,[sizes,1]);
yy = reshape(Come, [sizes, 1]);

b1 = xx\yy;
scatter(log10(xx),log10(yy))

vals = 0:1:3;
bs = vals*b1;
hold on
plot(vals,bs)
xlabel('log10 GITR Deposition Values','FontSize',axisLabelFont)
ylabel('log10 ERO Deposition Values','FontSize',axisLabelFont)
title('Least Squares Fit for GITR-ERO Comparison of Deposition','FontSize',titleFont)
set(gca,'FontSize',tickFont)
hold off
yCalc1 = b1*xx;
Rsq1 = 1 - sum((yy - yCalc1).^2)/sum((yy - mean(yy)).^2)
axis([0 3 0 3])