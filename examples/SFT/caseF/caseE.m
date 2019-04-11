clear all
close all

Ti0=linspace(20,50,100);
m=12;
mz=12;
Mi=0.1;
fcond = 1;
k0e=2000;
ne = 1e19;

cs = sqrt((2*1.602e-19*Ti0)/2/1.66e-27);
tau_s = m*Ti0.*(Ti0/2).^0.5./(6.8e4*(1+2/m)*ne/1e18*(4^2)*15);

FFf = 5*12*1.66e-27*0.1*cs./Ti0./tau_s./1.602e-19;
lambda_par = 1.602e-19*Ti0.*tau_s/mz/1.66e-27./(cs);
sinj=5;
sv=10;
delta_s = sv-sinj;
nz_np = exp(-Mi*delta_s./lambda_par.*(1-0.1*fcond/Mi));
figure(1)
plot(1./(Ti0.^2),FFf)
figure(2)
plot(1./(Ti0.^2),nz_np)
