clear all
close all
TD=50;
vTh = sqrt(50*1.602e-19/12/1.66e-27);
nP = 1000;
nT = 10000;
dt = 1e-8;

L=10;
sv = 1.2;
z0=0.1;
z1=0.2;
% odds = 1:2:nP;
z = (z1-z0)*rand(1,nP)+z0;
% z(odds) = z(odds) - 2*(z1-z0)+2*L-z0;
vz = zeros(1,nP);
vr = ones(1,nP);
% vz = normrnd(0,vTh/4,1,nP);
Dpar = 6.3357e3;
tau_s = 1.5756e-6;
tau_par = 1.8382e-6;
coeff_par = vTh*sqrt(2*dt/tau_par);
Dpar_distance = sqrt(Dpar*dt);
vE = dt*25*4*1.602e-19/12/1.66e-27;
dE = vE*dt;
vD = 1e4;
tau_s = 1.5756e-5;
v_therm = 2e8*dt;
nZ = 100;
gridZ = linspace(0,2*L,nZ);
dZ = gridZ(2) - gridZ(1);
spec = zeros(1,nZ);
inInd = 1:1:nP;
tic
for i=1:nT
    vtot = sqrt(vr.*vr + vz.*vz);
    outInds = find(z(inInd) < 0 | z(inInd)> 2*L);
    inInd(outInds) = [];
    %Diffusion
    r1 = normrnd(0,1,1,nP);
%     par = coeff_par./sqrt(abs(vTh*vTh - vz.*vz)).*r1.*abs(vz);
     par = coeff_par*r1;
%     %Drag
      drag = dt/tau_s*vz;
%     dragL = find(z < sv);
%     dragR = find(z > 2*L - sv);
%     vz(dragL) = vz(dragL) + (-vD - vz(dragL))/5e-6*dt;
%     vz(dragR) = vz(dragR) + (vD - vz(dragR))/5e-6*dt;
%     %thermal force
%     thermL = find(z < L);
%      thermR = find(z > L);
%      vz(thermL) = vz(thermL) + v_therm;
%      vz(thermR) = vz(thermR) - v_therm;
    diffz = abs(vz)./vtot.*par;
    diffz(find(vtot==0)) = 0;
    diffr = abs(vr)./vtot.*par;
    diffr(find(vtot==0)) = 0;
    vz = vz -vE-drag + diffz;
     vr = vr + diffr;
%     fast = find(abs(vz) > vTh);
%     vz(fast) = sign(vz(fast))*vTh;
    z(inInd) = z(inInd)+ vz(inInd)*dt;
    zinds = floor(z(inInd)./dZ)+1;
    zinds(find(zinds < 1)) = [];
    zinds(find(zinds > nZ)) = [];
    spec(zinds) = spec(zinds)+1;
%     if(mod(i,100) == 0)
%         figure(1)
%         semilogy(gridZ,spec)
%         a=1
%     end
end
toc
s = linspace(0,6);
figure(1)
semilogy(gridZ,1/7*spec)
hold on
plot(s,nP*exp(-2*s))

figure(3)
histogram(z)
