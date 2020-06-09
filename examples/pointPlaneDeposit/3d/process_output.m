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
Egrid = linspace(0,1000,200);
Agrid = linspace(0,90,30);
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
X = X(indices_to_plot,:);
Y = Y(indices_to_plot,:);
Z = Z(indices_to_plot,:);
grossDep = grossDep(indices_to_plot);
% patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),grossDep(surface),'FaceAlpha',.3)
patch(X',Y',Z',grossDep,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Deposited Impurity Mass (per face) [g]')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
axis tight
cb1 = colorbar
set(gca,'fontsize',16)
set(cb1,'fontsize',14);
% axis equal
    figure(3)
    edist2d = reshape(sum(edist,3),30,200);
h = pcolor(Egrid,Agrid,(edist2d))
h.EdgeColor = 'none';
colorbar
axis([5 1000 0 90])
 set(gca, 'XScale', 'log')

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

hold on
for i=1:100
plot3(x(:,i),y(:,i),z(:,i))
end