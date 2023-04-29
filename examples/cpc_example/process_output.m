clear all;
close all;
file = 'output/positions.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
weight = ncread(file,'weight');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
hitWall =ncread(file,'hitWall');
hit = find(hitWall);
not_hit = find(hitWall == 0);

E = 0.5*184*1.66e-27*(vx.^2 + vy.^2 + vz.^2)/1.602e-19;

figure(1)
scatter3(x(hit),y(hit),z(hit))
scatter3(x(not_hit),y(not_hit),z(not_hit))
title('Final W Particle Positions')
xlabel('X')
ylabel('Y')
zlabel('Z')

file = 'output/surface.nc';
grossDep = ncread(file,'grossDeposition');

grossEro = ncread(file,'grossErosion');

edist = ncread(file,'surfEDist');
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, Elen] = netcdf.inqDim(ncid,1);
[dimname, Alen] = netcdf.inqDim(ncid,2);
Egrid = linspace(0,1000,Elen);
Agrid = linspace(0,90,Alen);
sumWeightStrike = ncread(file,'sumWeightStrike');
aveSpyl = ncread(file,'aveSpyl');
spylCounts = ncread(file,'spylCounts');
sumParticlesStrike = ncread(file,'sumParticlesStrike');
totalW_WSputtYld = sum(aveSpyl(find(aveSpyl)))/sum(spylCounts(find(spylCounts)))
spylDist = ncread(file,'surfSputtDist');
rfylDist = ncread(file,'surfReflDist');
sumWeightStrike = ncread(file,'sumWeightStrike');
aveSpyl = ncread(file,'aveSpyl');
spylCounts = ncread(file,'spylCounts');
sumParticlesStrike = ncread(file,'sumParticlesStrike');


if (exist('x1') == 0)
    fid = fopen('input/gitrGeometryPointPlane3d.cfg');
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:19
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end

figure(2)
X = [transpose(x1),transpose(x2),transpose(x3)];
Y = [transpose(y1),transpose(y2),transpose(y3)];
Z = [transpose(z1),transpose(z2),transpose(z3)];

indices_to_plot = 1:1:length(X);
indices_to_plot(find(abs(y1+0.015)<0.0001 & abs(y2+0.0150)<0.001 & abs(y3+0.015)< 0.001)) = [];
% indices_to_plot(find(abs(z1-0.03)<0.0001 & abs(z2-0.03)<0.001 & abs(z3-0.03)< 0.001)) = [];
X = X(indices_to_plot,:);
Y = Y(indices_to_plot,:);
Z = Z(indices_to_plot,:);
grossDep = grossDep(indices_to_plot);
grossDep = 184/1e5*grossDep;
% patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),grossDep(surface),'FaceAlpha',.3)
patch(X',Y',Z',grossDep,'FaceAlpha',0.5,'EdgeAlpha', 1)%,impacts(surface)
title('Gross Deposited Impurity Mass [amu/s]')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis tight
cb1 = colorbar
set(gca,'ColorScale','log')
set(gca,'fontsize',16)
set(cb1,'fontsize',14);
% axis equal
    figure(3)
    edist2d = reshape(sum(edist,3),length(Agrid),length(Egrid));
h = pcolor(Egrid,Agrid,edist2d)
h.EdgeColor = 'none';
cb1 = colorbar
set(gca,'ColorScale','log')
axis([5 1000 0 90])
 set(gca, 'XScale', 'log')
 title({'Summed Surface W Ion','Energy-Angle Distribution'})
xlabel('Energy [eV]')
ylabel('Angle [deg]')
set(gca,'fontsize',16)
set(cb1,'fontsize',14);

file = 'output/history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
sizeArray = size(x);
nP = sizeArray(2);

% figure(2)
% hold on
% for i=1:10
% plot3(x(:,i),y(:,i),z(:,i),'r')
% end

specFile = 'output/spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');

figure(4)
slice1 = dens(:,:,:,end);
slice1 = reshape(slice1,length(gridR),length(gridY),length(gridZ));
slice1 = sum(slice1,2);
slice1 = reshape(slice1,length(gridR),length(gridZ));
dV = (gridR(2) - gridR(1))*(gridZ(2) - gridZ(1))*(gridY(2) - gridY(1))
slice1 = 1/1e5*1e-10/dV*slice1./length(gridY);
p = pcolor(gridR,gridZ,slice1')
p.EdgeColor = 'none';
hold on 
plot([-0.015, 0.015],[-0.00866025403784439,0.00866025403784439],'w','lineWidth',2)
axis equal
axis([-0.015 0.015 -0.01 0.02])
title({'W Impurity Density (All Charges)', 'Averaged in y-direction [m^{-3}]'})
xlabel('x [m]') % x-axis label
ylabel('z [m]') % y-axis label
set(gca,'fontsize',16)
colorbar
set(gca,'ColorScale','log')
% caxis([0 200000])