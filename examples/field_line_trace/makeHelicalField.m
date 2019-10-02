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
gridR = r + rC;
gridZ = z + zC;
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
ncid = netcdf.create(['./helixB.nc'],'NC_WRITE')
dimR = netcdf.defDim(ncid,'nR',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);
gridRvar = netcdf.defVar(ncid,'r','float',[dimR]);
gridZvar = netcdf.defVar(ncid,'z','float',[dimZ]);
xVar = netcdf.defVar(ncid,'bx','float',[dimR dimZ]);
yVar = netcdf.defVar(ncid,'by','float',[dimR dimZ]);
zVar = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);


netcdf.endDef(ncid);
netcdf.putVar(ncid, gridRvar, gridR);
netcdf.putVar(ncid, gridZvar, gridZ);
netcdf.putVar(ncid, xVar, br');
netcdf.putVar(ncid, yVar, by');
netcdf.putVar(ncid, zVar, bz');

netcdf.close(ncid);