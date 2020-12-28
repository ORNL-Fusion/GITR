%translate from 10,10,30 to 0,0,0
% clear all
% close all
% meshing2
%planes = planes./100;

planes(:,[1 2 4 5 7 8]) = (planes(:,[1 2 4 5 7 8]));
planes(:,[3 6 9]) = (planes(:,[3 6 9]));

materialPosZ = 0.0;
tol = 1e-6;
whereMat01 = find(abs(planes(:,3)) < tol);
whereMat02 = find(abs(planes(:,6)) < tol);
whereMat03 = find(abs(planes(:,9)) < tol);

whereMat = intersect(whereMat01,whereMat02);
whereMat = intersect(whereMat,whereMat03);
whereMat1 = find(planes(whereMat,1).^2 + planes(whereMat,2).^2 <= 0.049^2);
whereMat2 = find(planes(whereMat,4).^2 + planes(whereMat,5).^2 <= 0.049^2);
whereMat3 = find(planes(whereMat,7).^2 + planes(whereMat,8).^2 <= 0.049^2);
whereMat1 = intersect(whereMat1,whereMat2);
whereMat1 = intersect(whereMat1,whereMat3);
%%%Tower
whereTow01 = find(planes(:,3) > tol);
whereTow02 = find(planes(:,6) > tol);
whereTow03 = find(planes(:,9) > tol);

whereTow011 = find(planes(:,3) < 0.2);
whereTow021 = find(planes(:,6) < 0.2);
whereTow031 = find(planes(:,9) < 0.2);

whereTow = intersect(whereTow01,whereTow02);
whereTow11 = intersect(whereTow011,whereTow021);
whereTow11 = intersect(whereTow11,whereTow031);
whereTow = intersect(whereTow,whereTow03);
whereTow = intersect(whereTow,whereTow11);
whereTow1 = find(planes(whereTow,1).^2 + planes(whereTow,2).^2 <= 0.06^2);
whereTow2 = find(planes(whereTow,4).^2 + planes(whereTow,5).^2 <= 0.06^2);
whereTow3 = find(planes(whereTow,7).^2 + planes(whereTow,8).^2 <= 0.06^2);
whereTow1 = intersect(whereTow1,whereTow2);
whereTow1 = intersect(whereTow1,whereTow3);

materialZ = zeros(length(planes),1);
surfs = zeros(length(planes),1);

materialZ(whereMat(whereMat1)) = 74.0;
surfs(whereMat(whereMat1)) = 1.0;
surfs(whereTow(whereTow1)) = 1.0;
%Section to pick what to plot
r = sqrt(planes(:,1).^2 + planes(:,2).^2);
plotSet = find(r < .07 | (planes(:,1)>0 &  planes(:,2)<0))
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
figure(10)
patch(transpose(X),transpose(Y),transpose(Z),'green','FaceAlpha',.3,'EdgeColor',[0 0 0])%'none')
title('PISCES-A Simulated GITR Geometry')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('Boundary','Tungsten')
%patch(transpose([planes(:,1);planes(:,4);planes(:,7)]),transpose([planes(:,2);planes(:,5);planes(:,8)]),transpose([planes(:,3);planes(:,6);planes(:,9)]))
%patch(planes(:,[1 4 7])',planes(:,[2 5 8])',planes(:,[3 6 9])')
hold on
% X2 = X(whereMat(whereMat1),:);
% Y2 = Y(whereMat(whereMat1),:);
% Z2 = Z(whereMat(whereMat1),:);
X2 = [planes(whereMat(whereMat1),1),planes(whereMat(whereMat1),4),planes(whereMat(whereMat1),7)];
Y2 = [planes(whereMat(whereMat1),2),planes(whereMat(whereMat1),5),planes(whereMat(whereMat1),8)];
Z2 = [planes(whereMat(whereMat1),3),planes(whereMat(whereMat1),6),planes(whereMat(whereMat1),9)];
patch(transpose(X2),transpose(Y2),transpose(Z2),'blue','FaceAlpha',1,'LineStyle','none')

legend('Boundary','Tungsten')
axis equal
az = -65;
el = 45;
view(az, el);
figure(11)
X3 = [planes(whereTow(whereTow1),1),planes(whereTow(whereTow1),4),planes(whereTow(whereTow1),7)];
Y3 = [planes(whereTow(whereTow1),2),planes(whereTow(whereTow1),5),planes(whereTow(whereTow1),8)];
Z3 = [planes(whereTow(whereTow1),3),planes(whereTow(whereTow1),6),planes(whereTow(whereTow1),9)];
patch(transpose(X3),transpose(Y3),transpose(Z3),'blue','FaceAlpha',1,'LineStyle','none')
plane_coeffs