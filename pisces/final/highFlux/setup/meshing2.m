model = createpde;
importGeometry(model,'piscesA_scadFinal.stl');%'Block.stl'
pdegplot(model,'FaceLabels','on')
tic
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.005)%,'Hmax', 0.005)%, 'Hmin',0.004);%,'Hmin',0.5, 'Hmax',0.8 %0.001 used on pisces ,'Hmin',0.00001, 'Hmax',0.01
figure(1)
pdeplot3D(model,'FaceAlpha',0.5)

[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);
%volumePlanes = zeros(nTets*4,9);
%pointTally = zeros(nPoints,1);


%tess = delaunayn(rand(100,3));

%tess = sort(tess,2);
tess = transpose(t(1:4,:));%sort(t(1:4,:),1);
% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
       tess(:,[1 3 4]);tess(:,[2 3 4])];

% % find replicate faces
% faces = sortrows(faces);
% k = find(all(diff(faces)==0,1));
% 
% % delete the internal (shared) faces
% faces([k;k+1],:) = [];
% 
% surfacenodes = unique(faces(:));

faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);
%k = find((diff(faces)==0));
faces([k2;k2+1],:) = [];
% C = unique(faces, 'rows');
C = faces;
% figure(4)
% for i=1:length(C)
%     plot3(p(1,C(i,:)), p(2,C(i,:)), p(3,C(i,:)))
%     hold on
% end

planes = zeros(length(C), 9);

planes = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];
toc
% figure(3)
% planesR = sqrt((planes(:,1)-1).^2 + (planes(:,2)-1).^2);
% smallR = find(planesR < 0.06)
% scatter3([planes(smallR,1);planes(smallR,4); planes(smallR,7)],[planes(smallR,2);planes(smallR,5); planes(smallR,8)],[planes(smallR,3);planes(smallR,6); planes(smallR,9)])
%scatter3([planes(:,1);planes(:,4); planes(:,7)],[planes(:,2);planes(:,5); planes(:,8)],[planes(:,3);planes(:,6); planes(:,9)])