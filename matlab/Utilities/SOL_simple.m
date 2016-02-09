clear variables
D_perp = 1; %perpendicular diffusion coefficient (m/s)
L = 25; %Connection length
width = .05; %few cm - SOL width
Te = 25;
Ti = 25;
n0 = 1e19;
mi = 1.66e-27;
q = 1.602e-19;
cs = sqrt((Te+Ti)*q/mi);

lam_sol =sqrt(D_perp*L/cs); %characteristic SOL length

%Euler Method to solve for Mach numer along s

n_steps = 1e5;
grids = linspace(0,L,n_steps);
M = zeros(1,n_steps);
M(n_steps) = 0.999;
ds = L/(n_steps-1);

for z=(n_steps):-1:2
    M(z-1) = M(z) - ds*(1/L)*(1+M(z)^2)/(1-M(z)^2);
end

M = M - M(1);
M = M/M(n_steps);

plot(M)
 %1D relation of density to mach number for isothermal plasma
n = n0./(1+M.*M);
plot(n)

n_stepr = 1e2;
gridr = linspace(0,width,n_stepr);
dens = zeros(n_stepr,n_steps);

dens(1,:) = n;
for z=2:n_stepr
    dens(z,:) = n*exp(-z/n_stepr*width/lam_sol);
end
%surf(dens,'EdgeColor','none')


%%Temperature calc
T_t = 1;
cst = sqrt(2*T_t*q/mi);

q0 = 3.5*n0/2*cst*T_t*q;
k0e = 1.8e3; %W/mev^-7/2
k0i = 60;
Tse = (1+7*q0*grids./(2*k0e*T_t^(7/2))).^(2/7).*T_t;
Tsi = (1+7*q0*grids./(2*k0i*T_t^(7/2))).^(2/7).*T_t;
 plot(grids, Tsi)

 grad_Ts = q0/k0i*(1+7/2*q0/k0i*grids./Tsi).^(-5/7);
  plot(grids, grad_Ts)

  grad_Tr = -Ti/lam_sol*exp(-gridr./lam_sol);
    plot(gridr, grad_Tr)