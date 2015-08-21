

Tsr_e = (Tu^(7/2) - 7/2*q0e/k0e*r(1))^(2/7)*exp(-r(2)/lam_sol);
Tsr_i = (Tu^(7/2) - 7/2*q0i/k0i*r(1))^(2/7)*exp(-r(2)/lam_sol);

grad_Ts_e = q0e/k0e*(1+7/2*q0e/k0e*r(1)/Tsr_e)^(-5/7)*0*q;
grad_Ts_i = q0i/k0i*(1+7/2*q0i/k0i*r(1)/Tsr_i)^(-5/7)*q;

grad_Tr_e = -Tsr_e/lam_sol*q;
grad_Tr_i = -Tsr_i/lam_sol*q;
%{
Ns = 1e3;
Nr = 1e3;

grids = linspace(0,Lc,Ns);
gridr = linspace(0,width,Nr);

Tgride = zeros(Nr,Ns);
Tgridi = zeros(Nr,Ns);



Tgride(1,:) = (Tu^(7/2) - 7/2*q0/k0e*grids).^(2/7);
Tgride(2:end,:) = transpose(exp(-gridr(2:end)/lam_sol))*Tgride(1,:);



Tgridi(1,:) = (Tu^(7/2) - 7/2*q0/k0i*grids).^(2/7);
Tgridi(2:end,:) = transpose(exp(-gridr(2:end)/lam_sol))*Tgridi(1,:);

surf(Tgridi,'EdgeColor','none')
%}