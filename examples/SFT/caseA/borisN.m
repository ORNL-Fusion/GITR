clear all
close all

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

Ti = 10;
Te = 10;
n=1e19;

Energy = 10;
Z = 4;
amu = 12;
mz=amu;
mi =2;
vTh = sqrt(2*Ti*1.602e-19/mz/1.66e-27);
mu = mz/(mz+mi);
Beta_i = 3*(mu+5*sqrt(2)*Z^2*(1.1*mu^(5/2) - 0.35*mu^(3/2)) - 1)/(2.6-2*mu + 5.4*mu^2)

BMagPart = 1;
w = Z*Q*BMagPart/(amu*MI);
dt = 1e-8;% 0.34/w;
E = [0 0 0];
B = [0 0 1];
nP = 1e3;
vx = sqrt(2*Energy*Q/(MI*amu))*ones(nP,1);
vy = zeros(nP,1);
vz = zeros(nP,1);
v = [vx vy vz];
x = zeros(nP,1);
y = zeros(nP,1);
z = 0.1*rand(nP,1)+0.1;
r = [x y z];

nT = 100000;
% Constants used in Boris method Lorentz Integrator
q_prime = Z*Q/(amu*MI)*dt/2;
coeff = 2*q_prime/(1+(q_prime*BMagPart).^2);

vHist = zeros(nP,3,nT);
% vRelHist = zeros(nP,1,nT);
xHist = zeros(nP,3,nT);
vFlow = [0 0 -750]
dv_ITG = dt/amu/MI*Beta_i*1.38*1.602e-19;

%diagnostic
    netx0 = -0.03;
    netx1 = 0.03;
    nX = 100;
    netz0 = 0.0;
    netz1 = 2.0;
    nZ = 150;
    specBin = zeros(nX,nZ);
    specGridR = linspace(netx0,netx1,nX);
    specGridZ = linspace(netz0,netz1,nZ);
dR = specGridR(2) - specGridR(1);
    dZ = specGridZ(2) - specGridZ(1);
    specGridZ = linspace(netz0,netz1,nZ);
  tic  
% Boris Method Lorentz Integrator
for i=1:nT
%     vParticle = norm(v);
%     v = v0/vParticle*v;
    xind = floor((r(:,1) - netx0)/dR);
    yind = floor((r(:,3) - netz0)/dZ);
    xind(xind<=0) = [];
    yind(yind<=0) = [];
    xind(xind>nX) = [];
    yind(yind>nZ) = [];
    specBin(xind,yind) = specBin(xind,yind)+1;
    v_minus =v + q_prime*E;
    v = v_minus + q_prime*[v_minus(:,2)*B(3) - v_minus(:,3)*B(2), v_minus(:,3)*B(1) - v_minus(:,1)*B(3),v_minus(:,1)*B(2) - v_minus(:,2)*B(1)];
    
    v = v_minus + coeff*[v(:,2)*B(3) - v(:,3)*B(2), v(:,3)*B(1) - v(:,1)*B(3),v(:,1)*B(2) - v(:,2)*B(1)];
    
    v = v + q_prime*E;
   
    step = v*dt;
    
    r = r + step; 
    
    [tau_s nu_E] = drag(Ti,Te,n,Z,1,amu*MI,2*MI,v,vFlow);
    Dpar = vTh*vTh*tau_s;
    vrelative = v-vFlow;
    drag0 = -dt/tau_s*vrelative;
    %coeff_par = (2*floor(rand()+0.5) - 1)*sqrt(0.5*abs(nu_E)*dt*norm(vrelative)^2)*vrelative./norm(vrelative);
%     vrelative = v-vFlow;
%      vColl = drag*vrelative;
     drag0(r(3) > 1.2,3) = 0;
     vNew = v +  [0 0 dv_ITG];
     v = vNew;
     v(:,3) = v(:,3) + drag0(:,3);
     r(:,3) = r(:,3) + (2*floor(rand(nP,1)+0.5)-1).*sqrt(Dpar*dt);
     vHist(:,:,i) = v;
%     vRelHist(i) = norm(vrelative);
     xHist(:,:,i) = r;
end
toc
length(find(r(:,3)>0))
figure(1)
pcolor(specGridR,specGridZ,specBin')
figure(2)
plot(specGridZ,sum(specBin))
figure(3)
histogram(r(:,3))

figure(3)
plot3(xHist(1,:),xHist(2,:),xHist(3,:))

figure(2)
plot(xHist(1,:),xHist(3,:))
tau_s = drag(Ti,Te,n,Z,1,amu*MI,2*MI,v,vFlow)

figure(4)
plot(vRelHist)

figure(5)
plot(vHist(3,:))

E = 0.5*12*1.66e-27/1.602e-19*(vHist(1,:).^2+vHist(2,:).^2+vHist(3,:).^2);
figure(6)
plot(E)

function [tau_s nu_E] = drag(Ti,Te,n,z,zbackground,m,mbackground,v,vFlow)
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;
vRelative = v - vFlow;
v_norm = sqrt(vRelative(:,1).^2+vRelative(:,2).^2+vRelative(:,3).^2);

lam_d = sqrt(EPS0*Te/(n*Q));%only one q in order to convert to J
lam = 4*pi*n*lam_d^3;
gam = Q^4*z^2*(zbackground.^2)*log(lam)/(m*m*4*pi*EPS0*EPS0);
%gam2 = 0.238762895*z^2*zbackground^2*log(lam)/(m*m)*MI*MI;

a = mbackground./(2*Q*Ti); %q is just to convert units - no z needed

x = v_norm.^2.*a;
psi_prime = 2*sqrt(x./pi).*exp(-x);
psi_psiprime = erf(sqrt(x));
psi = psi_psiprime - psi_prime;
nu_0 = n*gam'*1./v_norm.^3;
nu_s = (1+m./mbackground)'.*psi.*nu_0;
tau_s = 1./nu_s;
nu_d = 2*(psi_psiprime - psi./(2*x)).*nu_0;
nu_par =  psi./x.*nu_0;
nu_E = 2*((m./mbackground)'.*psi - psi_prime).*nu_0;


end
