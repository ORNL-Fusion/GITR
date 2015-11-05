e0 = 8.85e-12;
me = 9.11e-31;
mi = 1.66e-27;
q = 1.602e-19;

B = [0 1.0 .14];

T=20;
n=1e19;
z=1;
zbackground=1;
m=mi;
mbackground=me;
v= [1e5 0 0]

flow_v = [-1e3 0 0]


lam_d = sqrt(e0*T/(n*zbackground^2*q));%only one q in order to convert to J
lam = 4*pi*n*lam_d^3;
gam = q^4*z^2*zbackground^2*log(lam)/(m*m*4*pi*e0*e0);
a = mbackground/(2*T*q);%q is just to convert units - no z needed
v_relative = v - flow_v;
x = norm(v_relative)^2*a;
psi_prime = 2*sqrt(x/pi)*exp(-x);
psi_psiprime = erf(sqrt(x));
psi = psi_psiprime - psi_prime;
nu_0 = gam*n/norm(v_relative)^3;
nu_s = (1+m/mbackground)*psi*nu_0;
nu_pitchangle = 2*(psi_psiprime - psi/(2*x))*nu_0;
nu_par = psi/x*nu_0;
nu_E = 2*(m/mbackground*psi - psi_prime)*nu_0;