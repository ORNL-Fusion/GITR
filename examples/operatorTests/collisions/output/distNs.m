close all
clear all
% load('nu.mat')
nP = 1e4;
% nT = 1e5;
m1 = 184*1.66e-27;
m2 = 2*1.66e-27;
nu_E = 2764;
k = 1.38e-23*11604; 
T=20;
v0 = 2000;

 meanSpeed = sqrt(2*k*T/(m2));
 meanSpeedImp = sqrt(2*k*T/(m1));
 v1 = collV(v0,meanSpeed,m1,m2);
% dv = abs(v0-v1);
nApprox = log(meanSpeedImp)/log(v1);% abs(meanSpeedImp - v0)/dv*(m2+m1)/m2
nBins = 100;


B = m2/(2*T*k);
vgrid = linspace(0,3*meanSpeed,nBins);
% fvb = nT*sqrt(B/pi)*exp(-B*vgrid.^2);
fvb = (B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
fvb = fvb/max(fvb);
fvbCDF = cumsum(fvb);
fvbCDF = fvbCDF./fvbCDF(end);
v = v0*ones(1,nP);
plot(vgrid,fvb)
run=1;
nT = 0;
vgrid2 = linspace(-20000,20000,nBins);
edges = linspace(-20000,20000,nBins+1);
B = m1/(2*T*k);
fvb2 = sqrt(B/pi)*exp(-B*vgrid2.^2);
fvb2 = fvb2/max(fvb2);
meas = 10;
tic
while meas > 1.15
    r1 = rand(1,nP);
    r2 = rand(1,nP);
    vb = interp1(fvbCDF,vgrid,r2,'pchip',0);
    phi= 2*rand(1,nP)-1;
    vb = phi.*vb;
    vf = collV(v,vb,m1,m2);
    v = vf;
    nT = nT+1;
    [N,edges] = histcounts(v,edges);
    N = N/max(N);
    meas = sum(sqrt((N - fvb2).^2));
end
toc
figure(6)
plot(vgrid2,fvb2)
hold on
plot(vgrid2,N)
B = m1/(2*T*k);
vx = linspace(-20000,20000);
fvb = nT*sqrt(B/pi)*exp(-B*vx.^2);
%fvb = nT*(B/pi)^(3/2)*4*pi*vx.*vx.*exp(-B*vx.^2);
vmag = sqrt(v.*v);
vy = v(randperm(length(v)));
vz = v(randperm(length(v)));
theta2 = 2*pi*rand(1,nP);
phi2 = pi*rand(1,nP);
% vx = v.*cos(theta2).*sin(phi2);
% vy = v.*sin(theta2).*sin(phi2);
% vz = v.*cos(phi2);
E = 0.5*m1*(v.*v +vy.*vy +vz.*vz)/1.602e-19;
vmag = sqrt(v.*v +vy.*vy +vz.*vz);
figure(2)
h1=histogram(v)
hold on
plot(vx,fvb*max(h1.Values)/max(fvb))

Egrid = linspace(0,200)*1.602e-19;
fe = 2*sqrt(Egrid./pi).*(1/(k*T))^(1.5).*exp(-Egrid./(k*T));
figure(3)
h2=histogram(E)
hold on
plot(Egrid/1.602e-19,fe*max(h2.Values)/max(fe))
nT
% 
% 
% v1 = 100;
% v2=100
% theta1 = -pi/4;
% theta2 = pi*5/4;
% m1=20;
% m2=20;
% phi = 0;
% coeff = (v1*cos(theta1-phi)*(m1-m2) + 2*m2*v2*cos(theta2-phi))/(m1+m2);
% 
% v1xPrime = coeff*cos(phi) + v1*sin(theta1-phi)*sin(phi+0.5*pi);
% v1yPrime = coeff*sin(phi) + v1*sin(theta1-phi)*cos(phi+0.5*pi);
function vf1 = collV(v1,v2,m1,m2)
vf1 = (m1-m2)./(m1+m2).*v1 + 2*m2/(m1+m2).*v2;
end