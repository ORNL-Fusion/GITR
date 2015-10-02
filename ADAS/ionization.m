%function z = ionization(z,dt,n,T,State,RateCoeff)

z=1;
T=750;
dt=1e-4;
n=1e11;

Nelem = size(Te(:,z));
diff = Te(:,z)-T;

for i=1:Nelem
    if (diff(i) >=0)
        i;
        break
    end
end

frac = (T-Te(i-1,z))/(Te(i,z)-Te(i-1,z));

Coeff = RateCoeff(i-1,z) + frac*RateCoeff(i,z)

Coeff2 = interp1(Te(:,z),RateCoeff(:,z),T)

if ( isnan(Coeff2) ) 
    error('Ionization interpolation out of range')
end

tion = 1/(Coeff*n)
P1 = exp(-dt/tion)


 if rand >= P1
    z = z+1
end   

