Eb = 9;
E = linspace(0,50,100);

f = E./(E+Eb).^3;
F = cumsum(f);
F = F/F(end);

%plot(E,F)

[val,ind] = min(abs(F-rand));

E(ind);

alph = linspace(-pi/2,pi/2,100);

fa = cos(alph);
Fa = cumsum(abs(fa));
Fa = Fa/Fa(end);

plot(alph, Fa)

[val2,ind2] = min(abs(Fa-rand));
[val3,ind3] = min(abs(Fa-rand));

Ka = [E(ind), alph(ind2)+pi/2, alph(ind3)+pi/2];