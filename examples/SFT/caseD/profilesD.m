close all
clear all
L = 10;
    nR = 10;
    r = linspace(-3,3,nR);
    
    
    
    nZ = 101;
    z = linspace(0.0,2*L,nZ);
    
    vz0 = -112.2;
    vz = zeros(1,nZ);
    vz(find(z<=1.2)) = vz0;
    vz2D = repmat(vz,nR,1);
    
    
    ncid = netcdf.create(['./profilesv112.nc'],'NC_WRITE')
    
    dimR = netcdf.defDim(ncid,'nR',nR);
    dimZ = netcdf.defDim(ncid,'nZ',nZ);
    
    gridRnc = netcdf.defVar(ncid,'gridR','double',dimR);
    gridZnc = netcdf.defVar(ncid,'gridZ','double',dimZ);
    Ez2Dnc = netcdf.defVar(ncid,'Ez','double',[dimR dimZ]);
    Ey2Dnc = netcdf.defVar(ncid,'Ey','double',[dimR dimZ]);
    Ex2Dnc = netcdf.defVar(ncid,'Ex','double',[dimR dimZ]);
    % Ni2Dnc = netcdf.defVar(ncid,'ni','double',[dimR dimZ]);
    % Te2Dnc = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
%     Ti2Dnc = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);
%     Te2Dnc = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
%     gradTi2Dnc = netcdf.defVar(ncid,'gradTiZ','double',[dimR dimZ]);
%     gradTe2Dnc = netcdf.defVar(ncid,'gradTeZ','double',[dimR dimZ]);
%     gradTi2DncR = netcdf.defVar(ncid,'gradTiR','double',[dimR dimZ]);
%     gradTe2DncR = netcdf.defVar(ncid,'gradTeR','double',[dimR dimZ]);
    
    %neVar = netcdf.defVar(ncid, 'Ne2', 'double',dimR);
    %teVar = netcdf.defVar(ncid, 'Te', 'double',dimR);
    netcdf.endDef(ncid);
    
    netcdf.putVar(ncid,gridRnc,r);
    netcdf.putVar(ncid,gridZnc,z);
    netcdf.putVar(ncid,Ez2Dnc,vz2D);
    netcdf.putVar(ncid,Ex2Dnc,0*vz2D);
    netcdf.putVar(ncid,Ey2Dnc,0*vz2D);
    
    
    %netcdf.putVar(ncid, neVar, Ne2);
%     netcdf.putVar(ncid, Ti2Dnc, Ti2D);
%     netcdf.putVar(ncid, Te2Dnc, Ti2D);
%     netcdf.putVar(ncid, gradTi2Dnc, gradTi);
%     netcdf.putVar(ncid, gradTe2Dnc, 0*gradTi);
%     netcdf.putVar(ncid, gradTi2DncR, 0*gradTi);
%     netcdf.putVar(ncid, gradTe2DncR, 0*gradTi);
    netcdf.close(ncid);
