
clear all
close all

Ti0=linspace(10,130,121);
m=12;
mz=m;
Z=4;
fcond = 1;
n0=1e19;
Mi=1;
k0e=2000;
mD=2;
mi=mD;
sinj=0.15;
sv=1.2;

phi_in = 1.73e23;
vTh = sqrt(Ti0*1.602e-19./m./1.66e-27);
tau_s = m*Ti0.*(Ti0/mD).^0.5./(6.8e4*(1+mD/m)*n0/1e18*(Z^2)*15);
vdiff = vTh.*vTh.*tau_s./sinj;
cs = sqrt((2*1.602e-19*Ti0)/mD/1.66e-27);
vB = -Mi*sqrt((2*1.602e-19*Ti0)/mD/1.66e-27);


mu = mz/(mz+mi);
Beta_i = 3*(mu+5*sqrt(2)*Z^2*(1.1*mu^(5/2) - 0.35*mu^(3/2)) - 1)/(2.6-2*mu + 5.4*mu^2)
s=sinj;
gamma=7;
mz=12;
mD=2
cs0=sqrt(2*Ti0*1.602e-19/mD/1.66e-27);
P=gamma*n0*cs0.*Ti0.*1.602e-19;
a = 7/2*fcond*P/k0e./Ti0.^(7/2);
dTids = Ti0.*(2*a./7./(a*s+1).^(5/7));
vTi = 1.602e-19/1.66e-27/m*Beta_i*tau_s.*dTids;
vpl = vdiff - vB - vTi;
figure(1)
plot(Ti0,phi_in./vpl./1e20,'lineWidth',2)
hold on
plot(Ti0,phi_in./vTh./1e20,'lineWidth',2)
title('Case B Peak Impurity Density')
xlabel('T_{i0} [eV]')
ylabel('n_p [10^{20} m^{-3}]')
axis([0 150 0 1])
set(gca,'fontsize',16)
legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}')
np=phi_in./vpl;




cs0=sqrt(2*Ti0*1.602e-19/mD/1.66e-27);
P=gamma*n0*cs0.*Ti0.*1.602e-19;
k0e=2000;
nS = 1000;
delta_s = (sv-sinj)/nS;
FFf=0*Ti0;
FiGf=0*Ti0;


for i=1:nS
    s = sinj+(i-1)*delta_s;
    Ts_s = Ti0.*(1+7/2*fcond*P*s/k0e./Ti0.^(7/2)).^(2/7);
    tau_ss = m*Ts_s.*(Ts_s/mD).^0.5./(6.8e4*(1+mD/m)*n0/1e18*(Z^2)*15);
    FFf = FFf+ mz*1.66e-27*(Mi*cs0)./Ts_s./1.602e-19./tau_ss.*delta_s;
    FiGf = FiGf + (Beta_i-1)*Ti0.*(2*a./7./(a*s+1).^(5/7))./Ts_s*delta_s;
end


nz_np = exp(-FFf+FiGf);
figure(2)
semilogy(1./(Ti0.^2)*1e4,nz_np.*np.*vTi./1.73e23)
% semilogy(1./(Ti0.^2)*1e4,exp(-40000./(Ti0.^2)))
axis([0 3.5 1e-4 1e-1])
nz_np.*np.*vTi./1.73e23;