cs = cs0*(Lc/(2*r(1)) - sqrt(Lc/(2*r(1)) - 1));
lam_sol =sqrt(D_perp*Lc/cs); %characteristic SOL length

ne_rs = n0*r(1)^2/Lc/(Lc/2 - sqrt((Lc/2)^2-r(1)^2))*exp(-(width -r(2))/lam_sol);
Tsr_e = (Tu^(7/2) - 7/2*q0e/k0e*r(1))^(2/7)*exp(-(width -r(2))/lam_sol);

if Tsr_e < C(1,1)
    frac = 1 - (C(1,1)  - Tsr_e)/C(1,1);
    rate = frac*C(1,2);
end
if Tsr_e > C(1,1)
diff = C(:,1) - Tsr_e;
ix = find(diff>0,1);
rate = C(ix-1,2) + (C(ix,1) - C(ix-1,1))/C(ix,1)*C(ix,2);
end

rate;
tion = 1/(rate*ne_rs);
P1 = exp(-dt/tion);

if rate > 0
 if rand >= P1
    ion = 1;
end   
end

