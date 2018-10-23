close all
clear all

m1 = 184*1.66e-27;
m2 = 2*1.66e-27;
nu_E = 2764;
k = 1.38e-23*11604; 
T=20;
v0 = 1000;

 meanSpeed = sqrt(2*k*T/(m2));
 meanSpeedImp = sqrt(2*k*T/(m1));
 v1 = collV(v0,meanSpeed,m1,m2);
% dv = abs(v0-v1);
nApprox = log(meanSpeedImp)/log(v1);% abs(meanSpeedImp - v0)/dv*(m2+m1)/m2
nT = 0;
dv = [];
% while v0<meanSpeedImp
for i=1:100
    v1 = collV(v0,meanSpeed,m1,m2);
    dv = [dv; v1-v0];
    v0 = v1;
    nT = nT+1;
end
nT
function vf1 = collV(v1,v2,m1,m2)
vf1 = (m1-m2)./(m1+m2).*v1 + 2*m2/(m1+m2).*v2;
end
