close all
clear all
load('nu.mat')
nP = 1e3;
nT = 1e3;
m1 = 184*1.66e-27;
m2 = 2*1.66e-27;
% nu_E = 2764;
dt = 1e-7;
v_b = 1000;
k = 1.38e-23*11604; 
T=20;
v0 = 3239;
% meanSpeed = sqrt(8*k*T/(m2*pi));
meanSpeed = sqrt(2*k*T/(m2));
% mean_dv = 2*m2/(m1+m2)*meanSpeed;
meanSpeedImp = sqrt(8*k*T/(m1*pi));
meanSpeedImp = sqrt(2*k*T/(m1));
v1 = collV(v0,0,m1,m2);
mean_dv = v0-v1;
collisions_per_acclimation = (abs(v0 - meanSpeedImp)/mean_dv);
percent = (nu_E*dt)*collisions_per_acclimation;
    if( percent >1)
        percent = 1;
    end
B = m2/(2*T*k);
vgrid = linspace(0,3*meanSpeed);
fvb = nT*sqrt(B/pi)*exp(-B*vgrid.^2);
fvb = nT*(B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
fvbCDF = cumsum(fvb);
fvbCDF = fvbCDF./fvbCDF(end);
v = v0*ones(1,nP);
vx = v;
vy = 0*vx;
vz = 0*vx;
plot(vgrid,fvb)
tic
% for i=1:nT
%     P1 = 1-exp(-nu_E*dt);
%     r1 = rand(1,nP);
%     y = find(r1 < P1);
%     r2 = rand(1,length(y));
%     vb = interp1(fvbCDF,vgrid,r2,'pchip',0);
%     phi= 2*rand(1,length(y))-1;
%     vb = phi.*vb;
%     vf = collV(v(y),vb,m1,m2);
%     v(y) = vf;
% end
% for i=1:nT
%     if(mod(i,nTapply)==1)
%     r2 = rand(1,nP);
%     vb = interp1(fvbCDF,vgrid,r2,'pchip',0);
%     phi= 2*rand(1,nP)-1;
%     vb = phi.*vb;
%     vf = collV(v,vb,m1,m2);
%     v = vf;
%     end
% end
E = 0.5*m1*(v.*v)/1.602e-19;
nu_E0 = interp1(Eparticle,nu_E,E,'pchip',-4.2116e5);
for i=1:nT
v1 = collV(v0,0,m1,m2);
mean_dv = v0-v1;
collisions_per_acclimation = (abs(v0 - meanSpeedImp)/mean_dv);
percent = abs(nu_E0(1)*dt*100);%abs(nu_E0./meanSpeed);%(nu_E0*dt)*collisions_per_acclimation;
    if( percent >1)
        percent = 1;
    end
    r1 = rand(1,nP);
    y = find(r1 < percent);
    r2 = rand(1,length(y));
    vb = interp1(fvbCDF,vgrid,r2,'pchip',0);
    phi= 2*rand(1,length(y))-1;
    vb = phi.*vb;
    vf = collV(v(y),vb,m1,m2);
    v(y) = vf;
end
toc

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