% fileID = fopen('Beers_Helicon_3D_100kW_HighDensity_Plasma-2.txt','r');
%fileID = fopen('Beers_Helicon_3D_100kW_LowDensity_Plasma.txt','r');
fileID = fopen(plasma_file,'r');
A = fscanf(fileID, '%f');
fclose(fileID);
A = reshape(A,17,[]);
Aorg=A;

extras = 102:102:102000;
A(:,extras) = [];
% A = A(:,1000:102000);
x = unique(A(1,:));
y = unique(A(2,:));
z = unique(A(3,:));

nX = length(x);
nY = length(y);
nZ = length(z);

bz = reshape(A(12,:),nY,nZ);
br = reshape(A(13,:),nY,nZ);
% bz = reshape(A(11,:),nY,nZ);

dens = reshape(A(14,:),nY,nZ);
electron_temp = reshape(A(15,:),nY,nZ);

ion_temp = reshape(A(16,:),nY,nZ);
slice = dens;
figure(1)
h = pcolor(y,z,dens')
h.EdgeColor = 'none';
colorbar

% figure
% h = pcolor(y,z,electron_temp')
% h.EdgeColor = 'none';
% colorbar


val = interpn(y,z,dens,0,0);
% set(gca, 'YDir', 'normal')
%  set(gca, 'XScale', 'log')
% title({'Sputtering Yield Kr on W','As a Function of Energy and Angle'})
% xlabel('E [eV]') % x-axis label
% ylabel('Angle [degrees]') % y-axis label
% set(gca,'fontsize',16)
% ax = h.Parent;   % Important
% set(ax, 'XTick', [ 0 50 100 200 400 600 1000])
% set(gca,'TickDir','out');
nX = 100;
nY = 110;
nZ = 200;
nR = 120;
r_grid = linspace(0.0,0.08,120);
x_grid = linspace(-0.08,0.08,100);
y_grid = linspace(-0.08,0.08,110);
z_grid = linspace(-0.8,0.8,200);

[rr_grid,zz_grid] = meshgrid(r_grid,z_grid);

% rr_grid = sqrt(xx_grid.^2 + yy_grid.^2);
% theta_grid = atan2(yy_grid,xx_grid);
br_interp = interpn(y,z,br,rr_grid,zz_grid);
bz_interp = interpn(y,z,bz,rr_grid,zz_grid);
bt_interp = 0*br_interp;


% bx_interp = br_interp.*cos(theta_grid);
% by_interp = br_interp.*sin(theta_grid);

te_interp = interpn(y,z,electron_temp,rr_grid,zz_grid);
ti_interp = interpn(y,z,ion_temp,rr_grid,zz_grid);
ne_interp = interpn(y,z,dens,rr_grid,zz_grid);

figure(11)
h = pcolor(r_grid,z_grid,bz_interp)
h.EdgeColor = 'none';
colorbar
title({'Proto-MPEX Magnetic Field (r) Profile'})
xlabel('r [m]')
ylabel('z [m]')
axis([0 max(y) min(z) max(z)])

br_interp(isnan(br_interp)) = 0;
bt_interp(isnan(bt_interp)) = 0;
bz_interp(isnan(bz_interp)) = 0;

te_interp(isnan(te_interp)) = 0;
ti_interp(isnan(ti_interp)) = 0;
ne_interp(isnan(ne_interp)) = 0;

%%

ncid = netcdf.create(['./profilesHelicon.nc'],'NC_WRITE')

dimR = netcdf.defDim(ncid,'nX',nR);
% dimY = netcdf.defDim(ncid,'nY',nY);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
% gridYnc = netcdf.defVar(ncid,'y','float',dimY);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
% % vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
% % vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
% % vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);
% % ernc = netcdf.defVar(ncid,'Er','float',[dimR dimZ]);
% % eznc = netcdf.defVar(ncid,'Ez','float',[dimR dimZ]);
% % etnc = netcdf.defVar(ncid,'Et','float',[dimR dimZ]);
% % gtirnc = netcdf.defVar(ncid,'gradTir','float',[dimR dimZ]);
% % gtiznc = netcdf.defVar(ncid,'gradTiz','float',[dimR dimZ]);
% % gtiync = netcdf.defVar(ncid,'gradTiy','float',[dimR dimZ]);
% % gternc = netcdf.defVar(ncid,'gradTer','float',[dimR dimZ]);
% % gteznc = netcdf.defVar(ncid,'gradTez','float',[dimR dimZ]);
% % gteync = netcdf.defVar(ncid,'gradTey','float',[dimR dimZ]);
% 
% 
% 
% %neVar = netcdf.defVar(ncid, 'Ne2', 'double',dimR);
% %teVar = netcdf.defVar(ncid, 'Te', 'double',dimR);
netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r_grid);
netcdf.putVar(ncid,gridZnc,z_grid);
netcdf.putVar(ncid,Ne2Dnc,ne_interp');
netcdf.putVar(ncid,Ni2Dnc,ne_interp');
netcdf.putVar(ncid,Te2Dnc,te_interp');
netcdf.putVar(ncid,Ti2Dnc,ti_interp');
% 
% % netcdf.putVar(ncid,vrnc,vrTotal');
% % netcdf.putVar(ncid,vtnc,vpTotal');
% % netcdf.putVar(ncid,vznc,vzTotal');
% 
netcdf.putVar(ncid,brnc,br_interp');
netcdf.putVar(ncid,btnc,bt_interp');
netcdf.putVar(ncid,bznc,bz_interp');
% 
% % netcdf.putVar(ncid,ernc,Epara');
% % netcdf.putVar(ncid,eznc,Eperp');
% % netcdf.putVar(ncid,etnc,(0*Epara)');
% % 
% % netcdf.putVar(ncid,gtirnc,gradTir');
% % netcdf.putVar(ncid,gtiznc,gradTiz');
% % netcdf.putVar(ncid,gtiync,gradTiy');
% % netcdf.putVar(ncid,gternc,gradTer');
% % netcdf.putVar(ncid,gteznc,gradTez');
% % netcdf.putVar(ncid,gteync,gradTey');
% 
% %netcdf.putVar(ncid, neVar, Ne2);
% %netcdf.putVar(ncid, teVar, Te);
netcdf.close(ncid);
