clear all
close all

nP=1e3;
nT=1e4;
v0 = 3e3;
flowV0 = 1e4;
T=20;
m1 = 184*1.66e-27;
m2 = 2*1.66e-27;
nu_E = 2764;
nu_d = 1000;
nu_par = 1000;
nu_friction = 1000;
dt = 1e-5;
k = 1.38e-23*11604; 

vx=zeros(1,nP);
vy=zeros(1,nP);
vz=v0*ones(1,nP);

flowVx=flowV0*ones(1,nP);
flowVy=zeros(1,nP);
flowVz=zeros(1,nP);

meanSpeed = sqrt(2*k*T/(m2));
meanSpeedImp = sqrt(2*k*T/(m1));

B = m2/(2*T*k);
vgrid = linspace(0,3*meanSpeed);
fvb = sqrt(B/pi)*exp(-B*vgrid.^2);
fvb = (B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
fvbCDF = cumsum(fvb);
fvbCDF = fvbCDF./fvbCDF(end);
figure(1)
plot(vgrid,fvb)
B_unit = [0,0,1];
tic
for i=1:nT
        if sum(isnan(vx))
        stop
    end
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
    n1 = normrnd(0,1,1,nP);
    n2 = normrnd(0,1,1,nP);
    xsi = rand(1,nP);
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
    vfx = collV(vx,vbx,m1,m2);
    vfy = collV(vy,vby,m1,m2);
    vfz = collV(vz,vbz,m1,m2);
    coeff_par = n1*sqrt(nu_par*dt);
    coeff_perp = abs(n2)*sqrt(nu_d*dt*0.5);
    cosXsi = cos(2.0*pi*xsi);
    sinXsi = sin(2.0*pi*xsi);
    vCollx = coeff_par.*d_parx + coeff_perp.*(perp_direction1(:,1)'.*cosXsi + perp_direction2(:,1)'.*sinXsi);
    vColly = coeff_par.*d_pary + coeff_perp.*(perp_direction1(:,2)'.*cosXsi + perp_direction2(:,2)'.*sinXsi);
    vCollz = coeff_par.*d_parz + coeff_perp.*(perp_direction1(:,3)'.*cosXsi + perp_direction2(:,3)'.*sinXsi);
    vCollNorm = sqrt(vCollx.*vCollx + vColly.*vColly + vCollz.*vCollz);
    vNew = sqrt(vfx.*vfx + vfy.*vfy + vfz.*vfz);
    thisSign = sign(vxRelative);
    vx = vxRelative + vRelativeNorm.*vCollx+ flowVx;
    vy = vyRelative + vRelativeNorm.*vColly;
    vz = vzRelative + vRelativeNorm.*vCollz;
     vPartNorm = sqrt(vx.^2 + vy.^2 + vz.^2);
     vx = vNew.*vx./vPartNorm - nu_friction*dt.*thisSign.*vRelativeNorm;
     vy = vNew.*vy./vPartNorm;
     vz = vNew.*vz./vPartNorm;
end
toc
B = m1/(2*T*k);
vgrid = linspace(0,20000);
fvb = nT*sqrt(B/pi)*exp(-B*vgrid.^2);
fvb = nT*(B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
vmag = sqrt((vx-0*flowV0).^2 +vy.*vy +vz.*vz);
figure(2)
h1=histogram(vmag)
hold on
plot(vgrid,fvb*max(h1.Values)/max(fvb))

vgrid = linspace(-20000,20000);
fvb = nT*sqrt(B/pi)*exp(-B*vgrid.^2);
figure(3)
h1=histogram(vx)
hold on
plot(vgrid,fvb*max(h1.Values)/max(fvb))
histogram(vy)
histogram(vz)

function vf1 = collV(v1,v2,m1,m2)
vf1 = (m1-m2)./(m1+m2).*v1 + 2*m2/(m1+m2).*v2;
end
