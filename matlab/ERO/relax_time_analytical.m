clear variables
T = 1e-6;
n_e = 1e19;
m(1) = 9.11e-31;
m(2) = 1.66e-27;

q = 1.602e-19;
cs = sqrt((2*T)*q/m(2));
e0 = 8.85e-12;

Z_test = 1;
M_test = m(2);
lam_d = sqrt(8.85e-12*T/(n_e*q));
lam = 4*pi*n_e*lam_d^3;
gam = q^4*Z_test^2*log(lam)/(M_test*M_test*4*pi*e0*e0);

v = 1e5*[1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
n_U = 1e4;
max = 3e4*cs;
U = linspace(0, max, n_U);
dU = U(2) - U(1);

a = m(1)/(2*T*q);
func = erf(U*sqrt(a))./U;
plot(U,func)

v_mag = norm(v);
v_mag = 1.39e4;
indx = round(v_mag/dU);

%dU_dt = n_e*gam*(func(indx+1) - func(indx-1))/(2*dU);

%tau_s = -v_mag/dU_dt

test_tau = v_mag^3/gam/(1+M_test/m(1))/n_e
