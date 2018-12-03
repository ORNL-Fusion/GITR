nn = 1;
nP = 1e4;
x = 5e3*ones(1,nP);
v = 1;
dt = 0.5e-6;

T=20;
m1 = 184*1.66e-27;
k = 1.38e-23*11604;% m2 kg s-2 eV-1
sigma = sqrt(k*T/m1); %meters per second
nu = 1e2;
D = sigma^2*nu;
nT = 2e4;
stepSize = sqrt(D*dt);

dt = dt/nn;
nT = nT*nn;
stepSize = stepSize/sqrt(nn);
tic
for i=1:nT
    r1 = rand(1,nP);
    plumin1 = 2*floor(r1 + 0.5)-1;
    x = x + stepSize*plumin1;
end
toc
figure(1)
h1 = histogram(x)
hold on
xgrid = linspace(-max(abs(x)),max(abs(x)));
% sigma = sigma;
mu = 0;
f = 1/sqrt(2*pi*sigma*sigma)*exp(-(xgrid-mu).^2/(2*sigma.^2));
plot(xgrid,f*max(h1.Values)/max(f))