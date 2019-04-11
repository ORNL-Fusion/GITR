clear all
close all

nP = 10000;
nT = 100000;
dt = 1e-8;

L=10;
sv = 1.2;
z0=0.1;
z1=0.2;
odds = 1:2:nP;
z = (z1-z0)*rand(1,nP)+z0;
z(odds) = z(odds) - 2*(z1-z0)+2*L-z0;
vz = zeros(1,nP);
Dpar = 3.6e3;

Dpar_distance = sqrt(Dpar*dt);
vD = 1e4;
v_therm = 2e8*dt;
nZ = 100;
gridZ = linspace(0,2*L,nZ);
dZ = gridZ(2) - gridZ(1);
spec = zeros(1,nZ);
inInd = 1:1:nP;
tic
for i=1:nT
    outInds = find(z(inInd) < 0 | z(inInd)> 2*L);
    inInd(outInds) = [];
    %Diffusion
    r1 = 2*floor(rand(1,nP)+0.5) - 1;
    z(inInd) = z(inInd)+ r1(inInd)*Dpar_distance;
    %Drag
    dragL = find(z < sv);
    dragR = find(z > 2*L - sv);
    vz(dragL) = vz(dragL) + (-vD - vz(dragL))/5e-6*dt;
    vz(dragR) = vz(dragR) + (vD - vz(dragR))/5e-6*dt;
    %thermal force
    thermL = find(z < L);
     thermR = find(z > L);
     vz(thermL) = vz(thermL) + v_therm;
     vz(thermR) = vz(thermR) - v_therm;
    
    z(inInd) = z(inInd)+ vz(inInd)*dt;
    zinds = floor(z(inInd)./dZ)+1;
    zinds(find(zinds < 1)) = [];
    zinds(find(zinds > nZ)) = [];
    spec(zinds) = spec(zinds)+1;
end
toc
figure(1)
semilogy(gridZ,spec)

figure(3)
histogram(z)
