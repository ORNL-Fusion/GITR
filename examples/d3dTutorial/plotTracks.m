clear all
close all
load('D3DrzZ2D2Tiles10Res.mat')
rGeom = r;
zGeom = z;
Zgeom = Z;


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


r = sqrt(x.^2 + y.^2);
figure(16)
hold on
for i=1:1000
plot(r(:,i),z(:,i))
end

plot(rGeom,zGeom)
hold on
scatter(rGeom,zGeom)
