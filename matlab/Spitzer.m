%function [Ts Td Tw] = Spitzer(z,zbackground,m,mbackground,v,n,T)
e0 = 8.85e-12;
me = 9.11e-31;
mi = 1.66e-27;
q = 1.602e-19;

B = [-0.2 1.0 0.3];

T=20;
n=1e19;
z=1;
zbackground=1;
m=12*mi;
mbackground=mi;
v= [0 0 1e4];

flow_v = [0 0 0];


lam_d = sqrt(e0*T/(n*zbackground^2*q));%only one q in order to convert to J
lam = 4*pi*n*lam_d^3;
gam = q^4*z^2*zbackground^2*log(lam)/(m*m*4*pi*e0*e0);
a = mbackground/(2*T*q);%q is just to convert units - no z needed


        B_unit = B/norm(B);
        g = v - flow_v;
        e3 = g/norm(g);
        
        s1 = dot(e3,B_unit);
        s2 = sqrt(1-s1^2);
        
        e1 = 1/s2*(s1*e3 - B_unit);
        e2 = -1/s2*cross(e3,B_unit);
%Spitzer slowing down time
v = v-flow_v;
x = sqrt(a)*norm(v);
G = (erf(x) - x*(2*exp(-x^2))/pi^(1/2))/(2*x^2);
G2 = erf(x) - x*(2/sqrt(pi)*exp(-x^2));
dU_dt = (1+m/mbackground)*gam*a*2*G*n;
tau_s = norm(v)*(1/((1+m/mbackground)*gam*a*2*G*n) )
tau_s2 = norm(v)^3/(gam*n*(1+m/mbackground)*G2);
%deflection time due to one species
dUU_dt_perp = 2*gam*n*(erf(x) - G)/norm(v);
tau_d = norm(v)^3*(1/(2*gam*n*(erf(x) - G)) )


%Energy exchange time

tau_e = norm(v)^3*(1/(4*2*gam*n*G));%%%%%Divide by a factor of 4 to get agreement with callen
W = 1/2*m*norm(v)^2;
dUU_dt = W/tau_e;

% function [e1 e2 e3] = direction(B, v, flow_v)
%         B_unit = B/norm(B);
%         g = v - flow_v;
%         e3 = g/norm(g);
%         
%         s1 = dot(e3,B_unit);
%         s2 = sqrt(1-s1^2);
%         
%         e1 = 1/s2*(s1*e3 - B_unit);
%         e2 = -1/s2*cross(e3,B_unit);
% end
A=1;
A1=1;
ts_slow = 11.7*A1^2*(T/8.62e-5)^(3/2)/((A+A1)*A1^(1/2)*n/1e6*z^2*zbackground^2*log(lam));% degrees in kelvin density in cm-3
ts_slow = 1.46e13*A1^2*(T)^(3/2)/((A+A1)*A1^(1/2)*n*z^2*zbackground^2*log(lam));% degrees in ev density in m-3
ts_slow = 1.46e13/1e18*A1^2*(T)^(3/2)/((A+A1)*A1^(1/2)*n/1e18*z^2*zbackground^2*log(lam));% degrees in ev density in m-3