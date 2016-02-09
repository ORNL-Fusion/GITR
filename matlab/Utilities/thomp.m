%function [energy_x, energy_y, energy_z] = thompson_vXYZ( bindingEnergy_eV, nP)

%Eb = bindingEnergy_eV;
nP=2;
Eb = 9;
maxE = 100

E = linspace(0,maxE,100);

f = E./(E+Eb).^3;
F = cumsum(f);
F = F/F(end);

plot(E,f)

[s1,s2,s3] = RandStream.create('mrg32k3a','NumStreams',3,'Seed','shuffle');

energy = interp1(F,E,rand(s1,1,nP))
%[val,ind] = min(abs(F-); % Energy



alph = linspace(-pi/2,pi/2,100);

fa = cos(alph);
Fa = cumsum(abs(fa));
Fa = Fa/Fa(end);

plot(alph, Fa)
stop
[val2,ind2] = min(abs(Fa-rand(s2,1,nP))); % Angle 1
[val3,ind3] = min(abs(Fa-rand(s3,1,nP))); % Angle 2

Ka = [E(ind), alph(ind2)+pi/2, alph(ind3)+pi/2];

energy_x = Ka(1);
energy_y = 1;
energy_z = 1;

%end