clear all
close all



model = createpde;
importGeometry(model,'helicon_6p256cm.stl');% Import STL file

figure(1)
pdegplot(model,'FaceLabels','on') %Plot stl

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.01);% Options Hmax and Hmin can be set, linear order can also be used

figure(2)
pdeplot3D(model,'FaceAlpha',0.5)

% Convert mesh to PET format (points, edges, triangles)
[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);

% Take 3D tetrahedral mesh and perform "skinning" operation
% Skinning operation takes the "outward facing triangles" from the
% mesh. This is used as the 3D surface mesh for GITR
tess = transpose(t(1:4,:));

% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
    tess(:,[1 3 4]);tess(:,[2 3 4])];


faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);

faces([k2;k2+1],:) = [];
% end of skinning operation

C = faces;

% planes is a variable which contains the xyz points of each triangle
% therefore it is a # triangles by 9 array
planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];
toc



plotSet = 1:1:length(planes);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

% Example of manual refinement
% this function splits the triangles in half N times, where N is the 
% 4th argument of the function
[X,Y,Z] = refineXYZ(X,Y,Z,1)
planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(3)
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','k')%'none')
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);

title({'Proto-MPEX Helicon Simulated GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')


nFaces = length(planes);
abcd = zeros(length(planes),4);
area = zeros(nFaces,1);
centroid = zeros(nFaces,3);
plane_norm = zeros(nFaces,1);
BCxBA = zeros(nFaces,1);
CAxCB = zeros(nFaces,1);


for i=1:length(planes)
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);
    
    AB = B-A;
    AC = C-A;
    BC = C-B;
    BA = A-B;
    CA = -AC;
    CB = -BC;
    
    norm1 = norm(AB);
    norm2 = norm(BC);
    norm3 = norm(AC);
    
    s = (norm1+norm2+norm3)/2;
    area(i) = sqrt(s*(s-norm1)*(s-norm2)*(s-norm3));
    normalVec = cross(AB,AC);
    
    d = -(dot(normalVec,A));
    
    abcd(i,:) = [normalVec,d];
    plane_norm(i) = norm(normalVec);
    
    BCxBA(i) = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i) = sign(dot(cross(CA,CB),normalVec));
    centroid(i,:) = [1/3*(planes(i,1)+planes(i,4)+planes(i,7)), ...
        1/3*(planes(i,2)+planes(i,5)+planes(i,8)), ...
        1/3*(planes(i,3)+planes(i,6)+planes(i,9))];
    rr = sqrt(centroid(1).^2 + centroid(2).^2);
end

% Material Z, surfaces, and surface normals have to be manually specified
% In this case, we use mathematical operations and assignment 
% on lines 138-143
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);


inDir = ones(1,length(planes));
figure(3)
hold on

for i=1:length(surfs)
    ind = i;
    normal = -abcd(ind,1:3)/plane_norm(ind);
    l_normal = 0.01;
    normal = l_normal*normal;

    dot_product = dot(normal,centroid(i,:));
    
    if(sign(dot_product) == 1)
        inDir(ind) = -1;
    end
    
    if(abs(normal(3)) > 0.005)
        surfs(ind) = 0;
        materialZ(ind) = 0;
    end
    
    normal = inDir(ind)*normal;
    
    %     quiver3(centroid(1),centroid(2),centroid(3),normal(1),normal(2),normal(3))
end

plotSet = find(surfs==0);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

% [X,Y,Z] = refineXYZ(X,Y,Z,6)
% planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(4)
patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',.3,'EdgeColor','k')%'none')
title('GITR non-surfaces: End-caps are absorbing boundaries')




format = '%.16e';

fileID = fopen('gitrGeometryPointPlane3d.cfg','w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,format,planes(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,1));
end
fprintf(fileID,' ] \n   y1 = [');
fprintf(fileID,format,planes(1,2));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,2));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,format,planes(1,3));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,3));
end
fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,format,planes(1,4));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,4));
end
fprintf(fileID,' ] \n   y2 = [');
fprintf(fileID,format,planes(1,5));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,5));
end
fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,format,planes(1,6));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,6));
end
fprintf(fileID,' ] \n   x3 = [');
fprintf(fileID,format,planes(1,7));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,7));
end
fprintf(fileID,' ] \n   y3 = [');
fprintf(fileID,format,planes(1,8));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,8));
end
fprintf(fileID,' ] \n   z3 = [');
fprintf(fileID,format,planes(1,9));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,9));
end
fprintf(fileID,' ] \n   a = [');
fprintf(fileID,format,abcd(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,1));
end

fprintf(fileID,' ] \n   b = [');
fprintf(fileID,format,abcd(1,2));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,2));
end

fprintf(fileID,' ] \n   c = [');
fprintf(fileID,format,abcd(1,3));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,3));
end
fprintf(fileID,' ] \n   d = [');
fprintf(fileID,format,abcd(1,4));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,4));
end
fprintf(fileID,' ] \n   plane_norm = [');
fprintf(fileID,format,plane_norm(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,plane_norm(i));
end
% fprintf(fileID,' ] \n   ABxAC = [');
% fprintf(fileID,format,ABxAC(1))
% for i=2:nFaces
% fprintf(fileID, ',')
% fprintf(fileID,format,ABxAC(i))
% end
fprintf(fileID,' ] \n   BCxBA = [');
fprintf(fileID,format,BCxBA(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,BCxBA(i));
end
fprintf(fileID,' ] \n   CAxCB = [');
fprintf(fileID,format,CAxCB(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,CAxCB(i));
end
fprintf(fileID,' ] \n   area = [');
fprintf(fileID,format,area(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,area(i));
end
fprintf(fileID,' ] \n   Z = [');
fprintf(fileID,'%f',materialZ(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%f',materialZ(i));
end
fprintf(fileID,' ] \n   surface = [');
fprintf(fileID,'%i',surfs(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%i',surfs(i));
end
fprintf(fileID,' ] \n   inDir = [');
fprintf(fileID,'%i',inDir(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%i',inDir(i));
end

fprintf(fileID,' ] \n');
fprintf(fileID,'periodic = 0;\n');
fprintf(fileID,'theta0 = 0.0;\n');
fprintf(fileID,'theta1 = 0.0\n');
fprintf(fileID,'periodic_bc_x0 = 0.0;\n');
fprintf(fileID,'periodic_bc_x1 = 0.0;\n');
fprintf(fileID,'periodic_bc_x = 0;}\n');
fclose(fileID)

function [Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n)

for j=1:n
    Xrefined = zeros(2*length(X),3);
    Yrefined = zeros(2*length(X),3);
    Zrefined = zeros(2*length(X),3);
    for i=1:length(X)
        A = [X(i,1) Y(i,1) Z(i,1)];
        B = [X(i,2) Y(i,2) Z(i,2)];
        C = [X(i,3) Y(i,3) Z(i,3)];
        
        AB = B-A;
        AC = C-A;
        BC = C-B;
        
        norms =[norm(AB) norm(BC) norm(AC)];
        [maxVal maxInd] = max(norms);
        if maxInd ==1
            midPtAB = A + 0.5*AB;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midPtAB(1);
            Xrefined(2*i-1,2) = midPtAB(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = C(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midPtAB(2);
            Yrefined(2*i-1,2) = midPtAB(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = C(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midPtAB(3);
            Zrefined(2*i-1,2) = midPtAB(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = C(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==2
            midptBC = B + 0.5*BC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = A(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = midptBC(1);
            Xrefined(2*i-1,3) = midptBC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = A(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = midptBC(2);
            Yrefined(2*i-1,3) = midptBC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = A(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = midptBC(3);
            Zrefined(2*i-1,3) = midptBC(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==3
            midptAC = A + 0.5*AC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midptAC(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = midptAC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midptAC(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = midptAC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midptAC(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = midptAC(3);
            Zrefined(2*i,3) = C(3);
        end
    end
    X = Xrefined;
    Y = Yrefined;
    Z = Zrefined;
end
end
