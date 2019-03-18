close all
clear all


TD0 = 10;

cs0 = sqrt((2*1.602e-19*TD0)/2/1.66e-27);


nR = 10;
r = linspace(-3,3,nR);



nZ = 101;
z = linspace(0.0,20.0,nZ);

flowVz = -.1*cs0;
vz = zeros(1,nZ);
vz(find(z<=10.85)) = flowVz;
vz2D = repmat(vz,nR,1);

gamma = 7;
n0=1e19;
fcond = 1;
k0e=2000;
P=gamma*n0*cs0.*TD0.*1.602e-19;
a = 7/2*fcond*P/k0e./TD0.^(7/2);
dTids = TD0.*(2*a./7./(a*z+1).^(5/7));

Ti = TD0*(1+7/2*fcond*P/k0e/TD0^(7/2)*z).^(2/7);
Ti2D = repmat(Ti,nR,1);
gradTi = repmat(dTids,nR,1);

ncid = netcdf.create(['./profiles.nc'],'NC_WRITE')

dimR = netcdf.defDim(ncid,'nR',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'gridR','double',dimR);
gridZnc = netcdf.defVar(ncid,'gridZ','double',dimZ);
vz2Dnc = netcdf.defVar(ncid,'vz','double',[dimR dimZ]);
vy2Dnc = netcdf.defVar(ncid,'vy','double',[dimR dimZ]);
vx2Dnc = netcdf.defVar(ncid,'vx','double',[dimR dimZ]);
% Ni2Dnc = netcdf.defVar(ncid,'ni','double',[dimR dimZ]);
% Te2Dnc = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);
gradTi2Dnc = netcdf.defVar(ncid,'gradTi','double',[dimR dimZ]);


%neVar = netcdf.defVar(ncid, 'Ne2', 'double',dimR);
%teVar = netcdf.defVar(ncid, 'Te', 'double',dimR);
netcdf.endDef(ncid);

netcdf.putVar(ncid,gridRnc,r);
netcdf.putVar(ncid,gridZnc,z);
netcdf.putVar(ncid,vz2Dnc,vz2D);
netcdf.putVar(ncid,vx2Dnc,0*vz2D);
netcdf.putVar(ncid,vy2Dnc,0*vz2D);


%netcdf.putVar(ncid, neVar, Ne2);
netcdf.putVar(ncid, Ti2Dnc, Ti2D);
netcdf.putVar(ncid, gradTi2Dnc, gradTi);
netcdf.close(ncid);