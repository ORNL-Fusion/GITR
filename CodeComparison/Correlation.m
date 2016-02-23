axisLabelFont = 20;
tickFont = 12;
titleFont = 36;

dA = CellProp{6}.data;
Come = NC*dA(1,1);
R = corrcoef(rot90(rot90(Tallys)),Come)

figure(11)
set(gcf,'Position',[300 300 900 600])
subplot(2,3,1)
surf(x,y,zeros(iNX,iNY),log10(rot90(rot90(Tallys))),'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([-1 3.5])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


subplot(2,3,4)
surf(x,y,zeros(iNX,iNY),log10(Come),'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([-1 3.5])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Deposition [log10(# particles)]','FontSize',titleFont)
set(gca,'FontSize',tickFont)


subplot(2,3,2)
surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_MQ)),'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([0 6])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Mean Charge [#]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,5)
surf(x,y,zeros(iNX,iNY),MQ,'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([0 6])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Mean Charge [#]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,3)
surf(x,y,zeros(iNX,iNY),rot90(rot90(gitr_ME)),'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([0 350])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('GITR Mean Impact Energy [eV]','FontSize',titleFont)
set(gca,'FontSize',tickFont)

subplot(2,3,6)
surf(x,y,zeros(iNX,iNY),ME,'EdgeColor','none')
colormap jet
colorbar
axis([-2 15 -5 8])
caxis([0 350])
xlabel('x axis [mm]','FontSize',axisLabelFont)
ylabel('y axis [mm]','FontSize',axisLabelFont)
zlabel('z axis [mm]','FontSize',axisLabelFont)
title('ERO Mean Impact Energy [eV]','FontSize',titleFont)
set(gca,'FontSize',tickFont)