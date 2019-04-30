% clear all
% close all
T=10;
n = 1e19;
dt=1e-7
tau1 = drag_tau(12,4,2,T,n)
vTh = sqrt(T*1.602e-19/12/1.66e-27);
0.5*12*1.66e-27*vTh*vTh/1.602e-19
D_par = vTh^2*tau1
tau_par = par_diff(12,T,4,2,T,n)
coeff_par = vTh*sqrt(2*dt/tau_par)
mfp = vTh*tau1
ds = vTh*tau1*dt
dfriction
TD = 100;
nD = 1e20;
mD = 2;
mi = 12;
Zi = 4;
fcond = 0.1;

tau = drag_tau(mi,Zi,mD,TD,nD)

vTh = sqrt(TD*1.602e-19/mi/1.66e-27);

D_par = vTh^2*tau
beta = Bi(mi,mD,Zi)
s=0.175;
dTids = gradTi(TD,mD,nD,fcond,s)

vTi = beta*tau*dTids*1.602e-19/mi/1.66e-27;
function tau_s = drag_tau(mi,Zi,mD,TD,nD)
nD = nD/1e18;
lnGam = 15;
tau_s = mi*TD*(TD/mD)^(1/2)/(6.8e4*(1+mD/mi)*nD*Zi^2*lnGam);
end

function tau_par = par_diff(mi,Ti,Zi,mD,TD,nD)
nD = nD/1e18;
lnGam = 15;
tau_par = mi*Ti*(TD/mD)^(1/2)/(6.8e4*nD*Zi^2*lnGam);
end

function Beta_i = Bi(mz,mi,Z)
mu = mz/(mz+mi);
Beta_i = 3*(mu+5*sqrt(2)*Z^2*(1.1*mu^(5/2) - 0.35*mu^(3/2)) - 1)/(2.6-2*mu + 5.4*mu^2);
end

function dTids = gradTi(Ti0,mD,n0,fcond,s)
gamma = 7;
k0e=2000;
cs0=sqrt(2*Ti0*1.602e-19/mD/1.66e-27);
P=gamma*n0*cs0.*Ti0.*1.602e-19;
a = 7/2*fcond*P/k0e./Ti0.^(7/2);
dTids = Ti0.*(2*a./7./(a*s+1).^(5/7));
end
