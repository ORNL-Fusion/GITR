clear all
close all

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;
dt = 1e-4;
zbackground = 1;
mbackground = 2*MI;
m2 =mbackground
nT =1e4;
nP = 1e3;
Te = 20;
Ti = 20;
T = Te;
n = 1e19;
z = 1;
k = 1.38e-23*11604; 
m = 184*MI;
m1 = m;
Eparticle = 10;% [ 0.1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 40 50 60 70 80 90 100];
thermal_speed = sqrt(2*Eparticle*Q/m);
cs = sqrt((Te+Ti)*Q*(mbackground.^-1));
v =  thermal_speed;
flowVelocity = 0;%[0 0 0];
v_relative = v;% - flowVelocity;
v_norm = v_relative;

     
lam_d = sqrt(EPS0*Te/(n*Q));%only one q in order to convert to J
lam = 4*pi*n*lam_d^3;
gam = Q^4*z^2*(zbackground.^2)*log(lam)/(m*m*4*pi*EPS0*EPS0);

dTheta2 = gam/2*n*dt/v_norm^3

v0 = 3e3;
flowV0 = 0;

vx=zeros(1,nP);
vy=zeros(1,nP);
vz=v0*ones(1,nP);
flowVx=flowV0*ones(1,nP);
flowVy=zeros(1,nP);
flowVz=zeros(1,nP);
B_unit = [1,0,0];
B = m2/(2*T*k);
meanSpeed = sqrt(2*k*T/(m2));
meanSpeedImp = sqrt(2*k*T/(m1));
vgrid = linspace(0,3*meanSpeed);
fvb = sqrt(B/pi)*exp(-B*vgrid.^2);
% fvb = (B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
fvbCDF = cumsum(fvb);
fvbCDF = fvbCDF./fvbCDF(end);
tic
for i=1:nT
    vPartNorm = sqrt(vx.^2 + vy.^2 + vz.^2);
    vxRelative = vx - flowVx;
    vyRelative = vy - flowVy;
    vzRelative = vz - flowVz;
    vRelativeNorm = sqrt(vxRelative.^2 + vyRelative.^2 + vzRelative.^2);
    d_parx = vxRelative./vRelativeNorm;
    d_pary = vyRelative./vRelativeNorm;
    d_parz = vzRelative./vRelativeNorm;
    s1 = d_parz;
    s2 = sqrt(1-s1.*s1);
    perp_direction1(:,1) = 1.0./s2.*(s1.*d_parx- B_unit(1));
    perp_direction1(:,2) = 1.0./s2.*(s1.*d_pary- B_unit(2));
    perp_direction1(:,3) = 1.0./s2.*(s1.*d_parz- B_unit(3));
    perp_direction2(:,1) = 1.0./s2.*(B_unit(2)*d_parz-B_unit(3)*d_pary);
    perp_direction2(:,2) = 1.0./s2.*(B_unit(3)*d_parx-B_unit(1)*d_parz);
    perp_direction2(:,3) = 1.0./s2.*(B_unit(1)*d_pary-B_unit(2)*d_parx);
    s20ind = find(s2==0 | abs(s1-1)<1.0e4);
    perp_direction1(s20ind,1) =  s1(s20ind);
    perp_direction1(s20ind,2) =  s2(s20ind);
    perp_direction1(s20ind,3) = 0;
    perp_direction2(s20ind,1) = d_parz(s20ind);
    perp_direction2(s20ind,2) = d_parx(s20ind);
    perp_direction2(s20ind,3) = d_pary(s20ind);
    s1(s20ind) = d_parx(s20ind)'.*perp_direction1(s20ind,1)+d_pary(s20ind)'.*perp_direction1(s20ind,2)+d_parz(s20ind)'.*perp_direction1(s20ind,3);
    s2(s20ind) = sqrt(1-s1(s20ind).*s1(s20ind));
    perp_direction1(s20ind,1) = -1.0./s2(s20ind)'.*(d_pary(s20ind)'.*perp_direction2(s20ind,3)-d_parz(s20ind)'.*perp_direction2(s20ind,2));
    perp_direction1(s20ind,2) = -1.0./s2(s20ind)'.*(d_parz(s20ind)'.*perp_direction2(s20ind,1)-d_parx(s20ind)'.*perp_direction2(s20ind,3));
    perp_direction1(s20ind,3) = -1.0./s2(s20ind)'.*(d_parx(s20ind)'.*perp_direction2(s20ind,2)-d_pary(s20ind)'.*perp_direction2(s20ind,1));

    r1 = rand(1,nP);
    r2 = rand(1,nP);
        r1x = rand(1,nP);
    r1y = rand(1,nP);
    r1z = rand(1,nP);
    r2 = rand(1,nP);
    vbx = interp1(fvbCDF,vgrid,r1x,'pchip',0);
    vby = interp1(fvbCDF,vgrid,r1y,'pchip',0);
    vbz = interp1(fvbCDF,vgrid,r1z,'pchip',0);
    phix= 2*rand(1,nP)-1;
    phiy= 2*rand(1,nP)-1;
    phiz= 2*rand(1,nP)-1;
    vbx = phix.*vbx;
    vby = phiy.*vby;
    vbz = phiz.*vbz;
    pm1= 2*rand(1,nP)-1;
    pm2= 2*rand(1,nP)-1;
    pm3= 2*rand(1,nP)-1;
    vdiffx = vx - vbx;
    vdiffy = vy - vby;
    vdiffz = vz - vbz;
    vdiffNorm = sqrt(vdiffx.*vdiffx + vdiffy.*vdiffy + vdiffz.*vdiffz);
    dTheta2 = gam/2*n*dt./vdiffNorm.^3;
    theta = sqrt(-2*dTheta2.*log(r1));
    phi= 2*pi*r2;
    du_par = vRelativeNorm.*cos(theta);
    du_perp1 = vRelativeNorm.*sin(phi).*sin(theta);
    du_perp2 = vRelativeNorm.*cos(phi).*sin(theta);
    mu = mbackground/(m+mbackground);
    dvx = mu*(du_par.*d_parx + du_perp1.*perp_direction1(:,1)' + du_perp2.*perp_direction2(:,1)');
    dvy = mu*(du_par.*d_pary + du_perp1.*perp_direction1(:,2)' + du_perp2.*perp_direction2(:,2)');
    dvz = mu*(du_par.*d_parz + du_perp1.*perp_direction1(:,3)' + du_perp2.*perp_direction2(:,3)');
    vx = vx + pm1.*dvx;
    vy = vy + pm2.*dvy;
    vz = vz + pm3.*dvz;
    a = 1;
end
toc
B = m1/(2*T*k);
vgrid = linspace(-20000,20000);
fvb = sqrt(B/pi)*exp(-B*vgrid.^2);
% fvb = (B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
figure(1)
h1 = histogram(vz)
hold on
plot(vgrid,fvb*max(h1.Values)/max(fvb))
figure(2)
hold on
histogram(vx)
histogram(vy)