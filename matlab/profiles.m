


%Background Plasma parameters
mb = mi;
Lc = 10; %Connection length
width = 0.05; %SOL width
%Define Fields
B0 = 1;
B = B0*[0 -.141 .99]; %Magnitude multiplied by a normalized vector
Tu = 25;
Tt = 1;
shtc_e = 7;
shtc_i = 3.5;
n0 = 1e19;
nt = n0/2;
cst = sqrt(2*Tt*q/mb);
k0e = 1.8e3; %W/mev^-7/2 heat conductivities for hydrogenic species
k0i = 60;
D_perp = 0.2;
Z = 1;
cs0 = sqrt((2*Tu)*q/mb);

q0e = shtc_e*nt*cst*Tt*q;
q0i = shtc_e*nt*cst*Tt*q;




%Test particle parameters
ion = 0; %ionization state 1 = true, 0=false
m = mi;
wc = q/m*B0;
Z_test = 1;
    mu = mb/(mb+m);
    alpha = Z_test^2*0.71;
    beta = -3*(1- mu - 5*sqrt(2)*Z_test^2)*(1.1*mu^(5/2) - 0.35*mu^(3/2))/(2.6 - 2*mu+ 5.4*mu^2);