ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 184;
Z = 1;

m = MI*amu;
q_m = Z*Q/m;
xp = [ 0 0 1e-3];
vp = [0 0 -3.2e3];
dt = 9e-8;
ld = 1e-5;
B = [0 0 -2];
BMagPart = norm(B);
E = [0 0 0];

surface_dz_dx = 1.73205;
surface_zIntercept = 0;
surfaceDirection = [surface_dz_dx 0 -1];
surfaceDirection_unit = surfaceDirection/norm(surfaceDirection);

dv_threshold = 0.001;

tic
v0 = vp;
r0 = xp;
perpDistanceToSurface = perpDist(r0,surface_dz_dx,surface_zIntercept);
while perpDistanceToSurface > 0

perpDistanceToSurface = perpDist(r0,surface_dz_dx,surface_zIntercept);
E = 45/(2*ld)*exp(-perpDistanceToSurface/(2*ld)) * surfaceDirection_unit;
k1r = dt*v0;
k1v = dt*q_m*(E + [v0(2)*B(3) - v0(3)*B(2), v0(3)*B(1) - v0(1)*B(3),v0(1)*B(2) - v0(2)*B(1)]);

perpDistanceToSurface = perpDist(r0+k1r/2,surface_dz_dx,surface_zIntercept);
E = 45/(2*ld)*exp(-perpDistanceToSurface/(2*ld)) * surfaceDirection_unit;
k2 = (v0+ k1v/2);
k2r = dt*k2;
k2v = dt*q_m*(E + [k2(2)*B(3) - k2(3)*B(2), k2(3)*B(1) - k2(1)*B(3),k2(1)*B(2) - k2(2)*B(1)]);


perpDistanceToSurface = perpDist(r0+k2r/2,surface_dz_dx,surface_zIntercept);
E = 45/(2*ld)*exp(-perpDistanceToSurface/(2*ld)) * surfaceDirection_unit;
k3 = (v0+ k2v/2);
k3r = dt*k3;
k3v = dt*q_m*(E + [k3(2)*B(3) - k3(3)*B(2), k3(3)*B(1) - k3(1)*B(3),k3(1)*B(2) - k3(2)*B(1)]);


perpDistanceToSurface = perpDist(r0+k3r,surface_dz_dx,surface_zIntercept);
E = 45/(2*ld)*exp(-perpDistanceToSurface/(2*ld)) * surfaceDirection_unit;
k4 = (v0+ k3v);
k4r = dt*k4;
k4v = dt*q_m*(E + [k4(2)*B(3) - k4(3)*B(2), k4(3)*B(1) - k4(1)*B(3),k4(1)*B(2) - k4(2)*B(1)]);

r1 = r0 + (k1r + 2*k2r + 2*k3r + k4r)/6;
v1 = v0 + (k1v + 2*k2v + 2*k3v + k4v)/6;
perpDistR1 = perpDist(r1,surface_dz_dx,surface_zIntercept);

if abs(norm(v1) - norm(v0))/norm(v0) > dv_threshold
                dt = dt/2;
else

v0 = v1;
r0 = r1;
end
end
toc
0.5*184*1.66e-27*norm(v1)^2/1.602e-19

                
