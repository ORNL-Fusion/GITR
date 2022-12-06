close all
clear all

file = strcat(pwd,'/history.nc');
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

for i=1:1

plot3(x(:,i),y(:,i),z(:,i),'k')
end
hold on
scatter3(x(1,1),y(1,1),z(1,1))
scatter3(x(end,1),y(end,1),z(end,1))
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Track','Begin','End')
axis equal
