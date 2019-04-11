clear all
close all
m=12;
mz=m;
Z=4;
mD = 2;
Mi = .1;
fcond = .1;
ne=1e20;
sinj=0.15;
sv=1.2;
delta_s = sv-sinj;

Ti0=linspace(10,130,100);

tau_s = m*Ti0.*(Ti0/mD).^0.5./(6.8e4*(1+mD/m)*ne/1e18*(Z^2)*15);

cs = sqrt((2*1.602e-19*Ti0)/mD/1.66e-27);

lambda_par = 1.602e-19*Ti0.*tau_s/mz/1.66e-27./(cs);

%from equation 6.80
nz_np = exp(-Mi*delta_s./lambda_par.*(1-0.1*fcond/Mi));

figure(2)
semilogy(1./(Ti0.^2)*1e4,nz_np,'lineWidth',2)
axis([0 3.5 1e-5 1])
xlabel('(T_t[eV])^{-2} x 10^{-4}')
ylabel('n_z(s)/n_p')