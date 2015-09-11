if r(2) <= lam_d

a = log(6)/lam_d;
phi0 = -1/2*Tu;

E_sh =  -phi0*a*exp(a*r(2));
else
    E_sh = 0;
end

Eext = [0 E_sh 0];