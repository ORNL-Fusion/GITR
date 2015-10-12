function [energy_x, energy_y, energy_z] = thompson_energyXYZ( bindingEnergy_eV, maxEnergy_eV, nP)

[s1,s2,s3] = RandStream.create('mrg32k3a','NumStreams',3,'Seed','shuffle');

% Energy

E = linspace(0,maxEnergy_eV,100);
f = E./(E+bindingEnergy_eV).^3;
F = cumsum(f);
F = F/F(end);

energy = interp1(F,E,rand(s1,1,nP));

% Angles
fudgeFac=1e-5
alph = linspace(-pi/2,pi/2-fudgeFac,100);
fa = cos(alph);
Fa = cumsum(abs(fa));
Fa = Fa/Fa(end);

theta = interp1(Fa,alph,rand(s2,1,nP));
phi = interp1(Fa,alph,rand(s3,1,nP));

energy_x = energy .* sin(theta) .* cos(phi);
energy_y = energy .* sin(theta) .* sin(phi);
energy_z = energy .* cos(theta);

end