%This data is for hydrogen on beryllium
E0 = 100;%incoming energy in eV - must be greater than threshold energy
alph = 10;%angle of incidence in degrees

%EB = ;%Surface binding energy
m1 = 1.66e-27;%projectile mass
m2 = 9*1.66e-27;%target mass
Q = 0.1;% Yield factor

Etf = 256;%Thomas-fermi energy
eps = E0/Etf;
Sn = 3.441*sqrt(eps)*log(eps+2.718)/(1+6.355*sqrt(eps)+eps*(6.882*sqrt(eps) - 1.708));

gam = 4*m1*m2/(m1+m2)^2;
Eth = 20;%EB/gam/(1-gam); %Threshold energy
del = Eth/E0;
g = (1-del^(2/3))*(1-del)^2;
YE0 = Q*Sn*g;

alph = alph*pi/180;
t = 1/cos(alph);
%f = ;
%sig = ;
YEa = YE0/t