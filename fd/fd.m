Te = 20;
Ti = 20;
q=1.602e-19;
e0 = 8.85e-12;
me = 9.11e-31;
mi = 1.66e-27;
d=0.1; %height of the simulation box
see = 0;%secondary electron emission coefficient for the material surface
ne = 2e19;

Dl =  sqrt(e0*Te/(ne*q));
gyr = 5e-3;


phi0 = -Te/(2*q)*log(2*pi*me/mi*(1+Ti/Te)/(1-see)^2)*q;

fd1 = 1- 