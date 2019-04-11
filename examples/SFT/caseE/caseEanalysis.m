close all
clear all

file = 'positions.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
hitWall =ncread(file,'hitWall');
hasLeaked =ncread(file,'hasLeaked');
nP=length(x);
hit = find(hitWall);
notHit = find(hitWall==0);
He_Wspyl = [4.07e-3 3.26e-3  3.02e-3];
He_Wroughness = [1.0 1.1 1.2];
E = 0.5*184*1.66e-27*(vx.^2 + vy.^2 + vz.^2)/1.602e-19;
Ehit = E(hit);
chargeHit = charge(hit);
r = sqrt(x.^2 + y.^2);
foundIndices = find(r<0.005);
figure(1)

scatter3(x,y,z)
hold on
%scatter3(x(hit),y(hit),z(hit))
title('Final W Particle Positions')
xlabel('X')
ylabel('Y')
zlabel('Z')

specFile = 'spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
% gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');
figure(2)
h=pcolor(gridR,gridZ,dens(:,:,5)')
h.EdgeColor = 'none';
figure(21)
semilogy(gridZ,sum(dens(:,:,5),1))
title('Impurity Density Profile')
xlabel('s [m]')
ylabel('n_p m^{-3}]')
set(gca,'fontsize',16)
maxDens = max(max(dens(:,:,5)));
area1 = find(gridZ< 2 | gridZ > 18);
area2 = find(gridZ> 2 & gridZ < 18);
totDens=sum(dens(:,:,5),1);

np = max(max(totDens(area1)));
npMiddle = max(max(totDens(area2)));

leakage = npMiddle/np;

file = 'history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');

leftSv = find(z>1 & z<1.2);
figure(101)
histogram(vz(leftSv))