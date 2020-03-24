close all
clear all

nP=1e6;
z0=0.1;
z1=0.2;
x0 = -0.02;
x1= 0.02;
x = (x1-x0)*rand(1,nP)+x0;
z = (z1-z0)*rand(1,nP)+z0;
y = 0.0*ones(1,nP);
T=10;
m=12;

vTh = sqrt(2*T*1.602e-19/m/1.66e-27);
k = 1.38e-23*11604; 
B = m*1.66e-27/(2*T*k);
vgrid = linspace(-3*vTh,3*vTh);
fv1 = sqrt(B/pi)*exp(-B*vgrid.^2);
fv1CDF = cumsum(fv1);
fv1CDF = fv1CDF./fv1CDF(end);
fvb = (B/pi)^(3/2)*4*pi*vgrid.*vgrid.*exp(-B*vgrid.^2);
fvbCDF = cumsum(fvb);
fvbCDF = fvbCDF./fvbCDF(end);
figure(9)
plot(vgrid,fvb)
vx = interp1(fv1CDF,vgrid,rand(1,nP),'pchip',0);
vy = interp1(fv1CDF,vgrid,rand(1,nP),'pchip',0);
vz = interp1(fv1CDF,vgrid,rand(1,nP),'pchip',0);
hold on
histogram(vx)
% phi  = rand(1,nP)*pi;
% theta = rand(1,nP)*2*pi;
% vz = vbx.*cos(phi);
% vx = vbx.*sin(phi).*cos(theta);
% vy = vbx.*sin(phi).*sin(theta);

vtot = sqrt(vx.^2 + vy.^2 + vz.^2);
histogram(vtot)
figure(10)
plot(vgrid,fv1)
figure(11)
histogram(vx)
hold on
histogram(vy)
histogram(vz)
vTh = sqrt(2*T*1.602e-19/m/1.66e-27);
0.5*12*1.66e-27*vTh*vTh/1.602e-19

% phi  = rand(1,nP)*pi;
% theta = rand(1,nP)*2*pi;
% vz = vTh*cos(phi);
% vx = vTh*sin(phi).*cos(theta);
% vy = vTh*sin(phi).*sin(theta);



% vx = vTh*ones(1,nP);
% vy = zeros(1,nP);
% vz = zeros(1,nP);

scatter(x,z)

ncid = netcdf.create(['./particleSourceD10Gauss.nc'],'NC_WRITE')
 
dimP = netcdf.defDim(ncid,'nP',nP);

xVar = netcdf.defVar(ncid,'x','double',[dimP]);
yVar = netcdf.defVar(ncid,'y','double',[dimP]);
zVar = netcdf.defVar(ncid,'z','double',[dimP]);
ExVar = netcdf.defVar(ncid,'vx','double',[dimP]);
EyVar = netcdf.defVar(ncid,'vy','double',[dimP]);
EzVar = netcdf.defVar(ncid,'vz','double',[dimP]);

netcdf.endDef(ncid);
 
netcdf.putVar(ncid, xVar, x);
netcdf.putVar(ncid, yVar, y);
netcdf.putVar(ncid, zVar, z);
netcdf.putVar(ncid, ExVar, vx);
netcdf.putVar(ncid, EyVar, vy);
netcdf.putVar(ncid, EzVar, vz);

netcdf.close(ncid);