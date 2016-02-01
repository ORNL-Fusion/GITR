function z = ionization(z,dt,n,T,State,RateCoeff)




z=1;
T=750;
dt=1e-9;
n=1e15;

Nelem = size(Te(:,z));
diff = Te(:,z)-T;

Coeff = interp1(Te(:,z),RateCoeff(:,z),T)

if ( isnan(Coeff) ) 
    error('Ionization interpolation out of range')
end

tion = 1/(Coeff*n)
P1 = 1-exp(-dt/tion)

r1=rand(s3)
 if r1 <= P1
    z = z+1
 end   


