function recombination(p,Bfield,Efield,xyz,dt,temp_eV,density,RecombRateCoeff,NTe,NNe)

z=1;
T=15;
dt=1e-4;
n=1e11;

            amu_mass = 1.66e-27;
            q = 1.602e-19;
            %T = 0.5*p.amu*amu_mass*(p.vx^2 +p.vy^2+ p.vz^2)/q;
            
            minT = min(DTEVD)
            minN = min(DDENSD)
            T=interpn(xyz.x,xyz.y,xyz.z,temp_eV(:,:,:,1),p.x,p.y,p.z)
            n=interpn(xyz.x,xyz.y,xyz.z,density(:,:,:,1),p.x,p.y,p.z)
            T = log10(T)
            n = log10(n)
           if (T > minT) && (n > minN)
            





Coeff = interp1(DTEVD,DDENSD,RecombRateCoeff(:,p.Z),T,n);

tion = 1/(Coeff*n)
P1 = exp(-dt/tion)


 if rand >= P1
    p.Z = p.Z-1
 end   
           end