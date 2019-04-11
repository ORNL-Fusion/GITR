close all
clear all
flowSpeed = [750 900 1000 1200 1400 1500];
L=10;
for i=1:length(flowSpeed)
    flowVz = -flowSpeed(i);
    
    TD0 = 10;
    dT_ds = 1.38;
    
    nR = 10;
    r = linspace(-300,300,nR);
    
    
    
    nZ = 101;
    z = linspace(0.0,2*L,nZ);
    
    vz = zeros(1,nZ);
    vz(find(z<=1.2)) = flowVz;
    vz(find(z>=2*L-1.2)) = -flowVz;
    vz2D = repmat(vz,nR,1);
    
    middle = ceil(nZ/2);
    Ti = zeros(1,nZ);
    Ti(1:middle) = 10+dT_ds*z(1:middle);
    Ti(middle+1:end) = 10+dT_ds*z(middle) + dT_ds*(L-z(middle+1:end));
    Ti2D = repmat(Ti,nR,1);
    
    gradTi = zeros(1,nZ);
    gradTi(1:middle) = dT_ds;
    gradTi(middle+1:end) = -dT_ds;
    gradTi2D = repmat(gradTi,nR,1);
    
    ncid = netcdf.create(['./profiles',num2str(flowSpeed(i)),'.nc'],'NC_WRITE')
    
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
    Te2Dnc = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
    gradTi2Dnc = netcdf.defVar(ncid,'gradTiZ','double',[dimR dimZ]);
    gradTe2Dnc = netcdf.defVar(ncid,'gradTeZ','double',[dimR dimZ]);
    gradTi2DncR = netcdf.defVar(ncid,'gradTiR','double',[dimR dimZ]);
    gradTe2DncR = netcdf.defVar(ncid,'gradTeR','double',[dimR dimZ]);
    
    
    
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
    netcdf.putVar(ncid, Te2Dnc, Ti2D);
    netcdf.putVar(ncid, gradTi2Dnc, gradTi2D);
    netcdf.putVar(ncid, gradTe2Dnc, 0*gradTi2D);
    netcdf.putVar(ncid, gradTi2DncR, 0*gradTi2D);
    netcdf.putVar(ncid, gradTe2DncR, 0*gradTi2D);
    netcdf.close(ncid);
end