exist_Te = exist('Te')
if exist_Te == 0
%background plasma parameters (maxwellian - one species)
T = 1;
n_e = 1e19;

%physical constants
mi = 1.66e-27;
me = 9.11e-31;
q = 1.602e-19;
e0 = 8.85e-12;

%Test particle parameters
m = mi;
Ka = [1, pi/4,pi/4]; 
v0 = sqrt(Ka(1)*2*q/m);
end
mi = 1.66e-27;
me = 9.11e-31;
n_e = 1e19;
T = Te;
e0 = 8.85e-12;
Z_test = 1;

%
lam_d = sqrt(e0*T/(n_e*q));
lam = 4*pi*n_e*lam_d^3;
gam = q^4*Z_test^2*log(lam)/(m*m*4*pi*e0*e0);
a = me/(2*T*q);


%Spitzer slowing down time
x = sqrt(a)*v0;
G = (erf(x) - x*(2*exp(-x^2))/pi^(1/2))/(2*x^2); 
tau_s = v0/((1+m/me)*gam*a*2*G*n_e)

%deflection time due to one species
tau_d = v0^3/(2*gam*n_e*(erf(x) - G))


%Energy exchange time
tau_e = v0^3/(4*2*gam*n_e*G)

