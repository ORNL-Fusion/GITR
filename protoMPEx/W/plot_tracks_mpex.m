clear all
close all
% M = csvread('iterGeom.csv');
% r = M(:,1);
% z = M(:,2);
% rGeom = r;
% zGeom = z;


file = 'output/history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
sizeArray = size(x);
nP = sizeArray(2);
hit = find(weight(end,:) < 1);

r = sqrt(x.^2 + y.^2);
figure(1)
hold on
% i=581
% plot(r(:,i),z(:,i))
%41242
for i=1:100:100000
 
% plot(r(:,i),z(:,i))
plot3(x(:,i),y(:,i),z(:,i))
end
xlabel('X')
ylabel('Y')
zlabel('Z')
% plot(rGeom,zGeom)
% hold on
% scatter(rGeom,zGeom)
axis equal
hold on
if (exist('x1') == 0)
    fid = fopen('input/gitrGeometry-ProtoMPEX-W.cfg');
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
end
Zsurface = Z;
surface = find(Zsurface);
nSurfaces = length(a);
%Section for finding a subset of planes to plot
r = sqrt(x1.^2 + y1.^2);
subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface)
title('Deposited Impurity Mass (per face) log scale [g]')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% axis([xMin xMax yMin yMax zMin zMax])
cb1 = colorbar
set(gca,'fontsize',16)
set(cb1,'fontsize',14);
axis equal
