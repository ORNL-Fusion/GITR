clear all;
close all;
base = 'output';

file = strcat(base,'/positions.nc');
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
weight = ncread(file,'weight');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
angle = ncread(file,'angle');
hitWall =ncread(file,'hitWall');
surfaceHit =ncread(file,'surfaceHit') + 1;
hit = find(hitWall);
not_hit = find(hitWall == 0);

E = 0.5*184*1.66e-27*(vx.^2 + vy.^2 + vz.^2)/1.602e-19;

figure
scatter3(x(hit),y(hit),z(hit))
hold on
scatter3(x(not_hit),y(not_hit),z(not_hit))
title('Final W Particle Positions')
xlabel('X')
ylabel('Y')
zlabel('Z')
%% Add charge resolution
row = histcounts(surfaceHit,0.5:1:12.5);

% Read in 3D surface mesh geometry file
if (exist('x1') == 0)
    fid = fopen(strcat(pwd,'/gitrGeometry.cfg'));

    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end
abcd = [a' b' c' d'];
surface = find(Zsurface);
nSurfaces = length(a);

subset = 1:length(x1);

figure
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),row,'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface)
title('Geometry')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar