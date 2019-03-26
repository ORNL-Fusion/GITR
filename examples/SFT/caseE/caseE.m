
clear all
close all

Ti0=linspace(10,200,191);
m=12;
mz=m;
Z=4;
fcond = .1;
n0=1e20;
Mi=.1;
k0e=2000;
mD=2;
mi=mD;
sinj=0.15;
sv=1.2;

phi_in = 1.73e23;
vTh = sqrt(Ti0*1.602e-19./m./1.66e-27);
tau_s = m*Ti0.*(Ti0/mD).^0.5./(6.8e4*(1+mD/m)*n0/1e18*(Z^2)*15);
vdiff = vTh.*vTh.*tau_s./.15;
cs = sqrt((2*1.602e-19*Ti0)/2/1.66e-27);
vB = -Mi*sqrt((2*1.602e-19*Ti0)/2/1.66e-27);


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




FFf = mz*1.66e-27*(Mi*cs0)./Ti0./1.602e-19./tau_s.*(sv-sinj);
Ts_sv = Ti0.*(1+7/2*fcond*P*sv/k0e./Ti0.^(7/2)).^(2/7);
Ts_sinj = Ti0.*(1+7/2*fcond*P*sinj/k0e./Ti0.^(7/2)).^(2/7);
FiGf = (Beta_i-1)*log(Ts_sv/Ts_sinj);
nz_np = exp(-FFf+FiGf);
figure(2)
semilogy(1./(Ti0.^2)*1e4,nz_np.*np.*vTi./1.73e23,'lineWidth',2)
% semilogy(1./(Ti0.^2)*1e4,exp(-40000./(Ti0.^2)))
axis([0 3.5 1e-5 1e-2])
xlabel('(T_0 [eV])^{-2} x 10^{-4}')
ylabel('\phi_{leak}/\phi_{in}')
title('Divertor Leakage as a Function of T^{-2}')
set(gca,'fontsize',16)
pbaspect([1 1.25 1])

figure(3)
M = csvread('../SFT-2/caseE-2-Table 1.csv',1,0)
Tgitr = M(:,1);
LeakGitr = M(:,6);
hold on 
scatter(1./(Tgitr.*Tgitr)*1e4,LeakGitr)
set(gca,'yscale','log')
axis([0 3.5 1e-4 2e-1])
xlabel('(T_0 [eV])^{-2} x 10^{-4}')
ylabel('\phi_{leak}/\phi_{in}')
title('Divertor Leakage as a Function of T^{-2}')
set(gca,'fontsize',16)
pbaspect([1 1.25 1])