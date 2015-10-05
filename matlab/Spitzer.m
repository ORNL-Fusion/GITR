%function [Ts Td Tw] = Spitzer(z,zbackground,m,mbackground,v,n,T)
e0 = 8.85e-12;
me = 9.11e-31;
mi = 1.66e-27;
q = 1.602e-19;

T=1e-5;
n=1e19;
z=1;
zbackground=1;
m=mi;
mbackground=me;
v= [1e5 0 0]


lam_d = sqrt(e0*T/(n*zbackground^2*q));%only one q in order to convert to J
lam = 4*pi*n*lam_d^3;
gam = q^4*Z^2*zbackground^2*log(lam)/(m*m*4*pi*e0*e0);
a = mbackground/(2*T*q);%q is just to convert units - no z needed



%Spitzer slowing down time
x = sqrt(a)*norm(v);
G = (erf(x) - x*(2*exp(-x^2))/pi^(1/2))/(2*x^2); 

tau_s = norm(v)*(1/((1+m/mbackground)*gam*a*2*G*n) )

%deflection time due to one species
tau_d = norm(v)^3*(1/(2*gam*n*(erf(x) - G)) )


%Energy exchange time
tau_e = norm(v)^3*(1/(4*2*gam*n*G))