%function z = ionization(z,dt,n,T,RecombState,RecombRateCoeff)

z=1;
T=15;
dt=1e-4;
n=1e11;

logn = log10(n);
logT = log10(T);

NTe = size(DTEVD);
NNe = size(DDENSD);
Tediff = DTEVD-logT;
Nediff = DDENSD-logn;


for i=1:NTe
    if (Tediff(i) >=0)
        i;
        break
    end
end

for j=1:NNe
    if (Nediff(j) >=0
        j;
        break
    end
end


lin1 = [RecombCoeff(i-1,j-1,z),]
Coeff = ;

tion = 1/(Coeff*n)
P1 = exp(-dt/tion)


 if rand >= P1
    z = z+1
end   