clear all
close all
rC = 1.6;
zC = 0.0;

r0 = 1;
r1 = 2.4;
z0 = -1.5
z1 = 1.5;

nR = 100;
nZ = 200;
br = zeros(nR,nZ);
by = ones(nR,nZ);
bz = zeros(nR,nZ);

r = linspace(-(r1-r0)/2,(r1-r0)/2,nR);
z = linspace(-(z1-z0)/2,(z1-z0)/2,nZ);

[r,z] = meshgrid(r,z);
angle = atan2(z,r);

br = 0.1*sin(angle);
bz = -0.1*cos(angle);
% br = 0.4*by;
% bz = -0.4*by;
% by = 2*by;
by = 0*br+1;
bmag = sqrt(br.^2 + bz.^2);
figure(1)
quiver(r,z,br,bz)

rP = 0.4;
zP = 0.1;
h = 0.01;
R = get_curvature(r,z,br,by,by,rP,zP)

function R = get_curvature(r,z,br,by,bz,rP,zP)
h = 0.01;
br_0 = interp2(r,z,br,rP,zP);
by_0 = interp2(r,z,by,rP,zP);
bz_0 = interp2(r,z,bz,rP,zP);
bmag = sqrt(br_0^2
br_1 = interp2(r,z,br,rP,zP);
R = br_0;
end