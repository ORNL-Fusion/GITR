clear all
close all
rC = 1.6;
zC = 0.0;

r0 = 1;
r1 = 2.4;
z0 = -1.5
z1 = 1.5;

nR = 1000;
nZ = 2000;
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
by = 0*br + 1;
bmag = sqrt(br.^2 + bz.^2);
% figure(1)
% quiver(r,z,br,bz)

rP = 0.1;
yP = 0.0
zP = 0.1;
h = 0.01;
R = get_curvature(r,z,br,by,bz,rP,yP,zP)
R=0.1412;
sqrt(rP^2+zP^2)
s = 0.02;
theta = linspace(0,pi,1000);
rr = sqrt(rP^2 + zP^2);
rr0 = sqrt(rP^2 + zP^2)-s;
rr1 = sqrt(rP^2 + zP^2)+s;

f_theta = (2*R*theta -s*sin(theta))/(2*pi*R);
f_theta = (2*R*s*theta -s*s*sin(theta)+R*R - (R-s)^2)/(2*pi*R*s);
f_theta = pi*(2*R*s-s*s*cos(theta));
% f_theta = linspace((rr^2-rr0^2),rr1^2-rr^2,length(theta))
% f_theta = 0*theta +1;
% plot(theta,2*pi*R*theta)
% hold on
% plot(theta,s*sin(theta))
plot(theta,f_theta)
cdf = cumsum(f_theta);
cdf = cdf./cdf(end);
plot(theta,cdf)

drand = rand(1,1e6);
theta_new = interp1(cdf,theta,drand);
nInside = length(find(theta_new < 0.5*pi));
nOutside = length(find(theta_new > 0.5*pi));



dens0 = nInside/(rr^2-rr0^2);
dens1 = nOutside/(rr1^2-rr^2);
function R = get_curvature(r,z,br,by,bz,rP,yP,zP)
h = 0.01;
br_0 = interp2(r,z,br,rP,zP);
by_0 = interp2(r,z,by,rP,zP);
bz_0 = interp2(r,z,bz,rP,zP);
bmag = sqrt(br_0^2 + by_0^2 + bz_0^2);
r_plus = rP+br_0/bmag*h;
y_plus = yP+by_0/bmag*h;
z_plus = zP+bz_0/bmag*h;
r_minus = rP-br_0/bmag*h;
y_minus = yP-by_0/bmag*h;
z_minus = zP-bz_0/bmag*h;

br_plus = interp2(r,z,br,r_plus,z_plus);
by_plus = interp2(r,z,by,r_plus,z_plus);
bz_plus = interp2(r,z,bz,r_plus,z_plus);
bmag_plus = sqrt(br_plus^2 + by_plus^2 + bz_plus^2);

br_minus = interp2(r,z,br,r_minus,z_minus);
by_minus = interp2(r,z,by,r_minus,z_minus);
bz_minus = interp2(r,z,bz,r_minus,z_minus);
bmag_minus = sqrt(br_minus^2 + by_minus^2 + bz_minus^2);

br_deriv2 = (br_plus/bmag_plus - 2*br_0/bmag + br_minus/bmag_minus)/(h^2); 
by_deriv2 = (by_plus/bmag_plus - 2*by_0/bmag + by_minus/bmag_minus)/(h^2); 
bz_deriv2 = (bz_plus/bmag_plus - 2*bz_0/bmag + bz_minus/bmag_minus)/(h^2);

br_deriv1 = (br_plus - br_minus)/(2*h);
by_deriv1 = (by_plus - by_minus)/(2*h); 
bz_deriv1 = (bz_plus - bz_minus)/(2*h);

cross1 = cross([br_deriv1, by_deriv1, bz_deriv1],[br_deriv2, by_deriv2, bz_deriv2]);
denom = norm([br_deriv1, by_deriv1, bz_deriv1]);

R = 1/(norm(cross1)/denom^3);
end