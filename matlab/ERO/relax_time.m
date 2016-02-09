clear variables
n_v = 100;
T = 25;
n_e = 1e19;
mi = 1.66e-27;
q = 1.602e-19;
cs = sqrt((2*T)*q/mi);
max = 3*cs;
v_grid = linspace(-max, max, n_v);
dv = v_grid(2) - v_grid(1);
dist = zeros(1,n_v);


    dist = exp(-mi*v_grid.*v_grid/(2*T*q));
    C = (mi/(2*pi*T*q))^(3/2);

plot(dist)

integr = (sum(dist)*dv)^3*C

Z_test = 1;
M_test = mi;
lam_d = sqrt(8.85e-12*T/(n_e*q));
lam = 4*pi*n_e*lam_d^3;
gam = 4*pi*q^4*Z_test^2*log(lam)/(M_test*M_test);
U = [0 0 0];
dU = [10 0 0];
U_prime = U+dU;

integrand = zeros(n_v,3);
integrand_p = zeros(n_v,3);
Nt = 1000;
int_hist = zeros(Nt,2);
for z=1:Nt

integrand(:,1) = dist./abs(-v_grid+U(1));
integrand(:,2) = dist./abs(-v_grid+U(2));
integrand(:,3) = dist./abs(-v_grid+U(3));

%integrand_p(:,1) = dist./abs(-v_grid+U_prime(1));
%integrand_p(:,2) = dist./abs(-v_grid+U_prime(2));
%integrand_p(:,3) = dist./abs(-v_grid+U_prime(3));

integr2 = sum(integrand)*dv;
integr2 = prod(integr2)*C;

%integr2_p = sum(integrand_p)*dv;
%integr2_p = prod(integr2_p)*C;
int_hist(z,1) = U(1);
U(1) = U(1) + dU(1);
int_hist(z,2) = integr2;
end
%dU_dt = n_e*gam*(M_test+mi)/mi*(integr2_p - integr2)/dU(1);

%ts = -U/dU_dt
plot(int_hist(:,1),int_hist(:,2))