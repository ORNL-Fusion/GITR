ne_rs = n0*r(1)^2/Lc/(Lc/2 - sqrt((Lc/2)^2-r(1)^2))*exp(-(width - r(2))/lam_sol);

lam_d = sqrt(e0*Tsr_e/(ne_rs*q));
lam = 4*pi*ne_rs*lam_d^3;
gam = q^4*Z_test^2*log(lam)/(m*m*4*pi*e0*e0);
a = me/(2*Tsr_e*q);
A = mb/(2*Tsr_i*q);


%Spitzer slowing down time
x = sqrt(a)*norm(v);
x2 = sqrt(A)*norm(v);
G = (erf(x) - x*(2*exp(-x^2))/pi^(1/2))/(2*x^2); 
G2 = (erf(x2) - x2*(2*exp(-x2^2))/pi^(1/2))/(2*x2^2); 
tau_s = norm(v)*(1/((1+m/me)*gam*a*2*G*ne_rs) + 1/((1+m/mb)*gam*A*2*G2*ne_rs));

%deflection time due to one species
tau_d = norm(v)^3*(1/(2*gam*ne_rs*(erf(x) - G)) + 1/(2*gam*ne_rs*(erf(x2) - G2)));


%Energy exchange time
tau_e = norm(v)^3*(1/(4*2*gam*ne_rs*G)+ 1/(4*2*gam*ne_rs*G2));

if tau_d < 0
    tau_d = 1e10;
    
end
if tau_e < 0
    tau_e = 1e10;
end
if tau_s < 0
    tau_s = 1e10;
    
end