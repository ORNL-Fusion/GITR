close all
clear all

filename = 'W Target Assembly-1.stl';
max_side_length = 0.5;
geometryFilename = 'gitrGeometry.cfg';

model = createpde;
importGeometry(model,filename);
pdegplot(model,'FaceLabels','on')
tic
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',max_side_length);
figure(1)
pdeplot3D(model,'FaceAlpha',0.5)

[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);

tess = transpose(t(1:4,:));
% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
       tess(:,[1 3 4]);tess(:,[2 3 4])];

% % find replicate faces
faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);

faces([k2;k2+1],:) = [];

C = faces;

planes = zeros(length(C), 9);

planes = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];
toc

materialZ = ones(1,length(planes));
materialZ(700:end)=0;
surfs = ones(1,length(planes));
surfs(700:end) = 0;

plotSet = find(surfs);
notSet = find(surfs == 0);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
figure(10)
patch(transpose(X),transpose(Y),transpose(Z),'green','FaceAlpha',.3,'EdgeColor',[0 0 0])%'none')
title('PISCES-A Simulated GITR Geometry')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('Boundary','Tungsten')

hold on

X2 = [planes(notSet,1),planes(notSet,4),planes(notSet,7)];
Y2 = [planes(notSet,2),planes(notSet,5),planes(notSet,8)];
Z2 = [planes(notSet,3),planes(notSet,6),planes(notSet,9)];
patch(transpose(X2),transpose(Y2),transpose(Z2),'blue','FaceAlpha',1,'LineStyle','none')

legend('Boundary','Tungsten')
axis equal
az = -65;
el = 45;
view(az, el);

GITR_write_3d_geometry_file