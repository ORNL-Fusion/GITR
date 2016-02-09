v = linspace(-10,10,10000);
dv = v(2) - v(1);
U = 1;
a = 1e4;
y = exp(-a*(v).*(v))./(abs(U-v)).^3.*(U-v);

plot(v,y)

int = sum(y)*dv

answer = 1/U*(pi/a)^(1/2)*erf(U*sqrt(a))