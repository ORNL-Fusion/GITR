
E_par = Tsr_e*(Lc/(2*r(1)*sqrt((Lc/2)^2 - r(1)^2)) - 1/r(1));

fs = log(Lc+sqrt(Lc^2-4*r(1)^2)) - log(2*Lc);
E_perp = -grad_Tr_e*fs;