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
T=100;
m=12;
vTh = sqrt(2*T*1.602e-19/m/1.66e-27);
0.5*12*1.66e-27*vTh*vTh/1.602e-19

phi  = rand(1,nP)*pi;
theta = rand(1,nP)*2*pi;
vz = vTh*cos(phi);
vx = vTh*sin(phi).*cos(theta);
vy = vTh*sin(phi).*sin(theta);



% vx = vTh*ones(1,nP);
% vy = zeros(1,nP);
% vz = zeros(1,nP);

scatter(x,z)

ncid = netcdf.create(['./particleSourceD100.nc'],'NC_WRITE')
 
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