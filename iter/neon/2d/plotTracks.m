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


r = sqrt(x.^2 + y.^2);
figure(1)
hold on
% i=581
% plot(r(:,i),z(:,i))
for i=1:100
plot(r(:,i),z(:,i))
% plot3(x(:,i),y(:,i),z(:,i))
end
xlabel('X')
ylabel('Y')
zlabel('Z')
% plot(rGeom,zGeom)
% hold on
% scatter(rGeom,zGeom)
axis equal
readIterGeom