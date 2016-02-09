clear variables
D_perp = 1; %perpendicular diffusion coefficient (m/s)
L = 25; %Connection length
width = .05; %few cm - SOL width
Te = 25;
Ti = 25;
n0 = 1e19;
mi = 1.66e-27;
m = mi;
Z_test = 1;
q = 1.602e-19;
cs = sqrt((Te+Ti)*q/mi);

s = 12;

lam_sol =sqrt(D_perp*L/cs); %characteristic SOL length

n_steps = 1e3;
grids = linspace(0,L,n_steps);

n_stepr = 1e2;
gridr = linspace(0,width,n_stepr);
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

  grad_Tr = -Tsi(round(s/L*n_steps))/lam_sol*exp(-gridr./lam_sol);
    plot(gridr, grad_Tr)
    
    mu = mi/(mi+m);
    alpha = Z_test^2*0.71;
    beta = -3*(1- mu - 5*sqrt(2)*Z_test^2)*(1.1*mu^(5/2) - 0.35*mu^(3/2))/(2.6 - 2*mu+ 5.4*mu^2);
    
   % F_therms = alpha*grad