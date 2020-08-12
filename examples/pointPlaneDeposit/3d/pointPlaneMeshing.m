clear all
close all
model = createpde;
importGeometry(model,'pointPlane.stl');% Import STL file
figure(2)
pdegplot(model,'FaceLabels','on') %Plot stl 

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.03);% Options Hmax and Hmin can be set, linear order can also be used
figure(1)
pdeplot3D(model,'FaceAlpha',0.5)

[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);

tess = transpose(t(1:4,:));%sort(t(1:4,:),1);
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

C = faces;

planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];
toc






% tol = 1e-6;
% verticalInds1 = find(abs(planes(:,1) - planes(:,4)) < tol & abs(planes(:,1) - planes(:,7)) < tol);
% verticalInds2 = find(abs(planes(:,2) - planes(:,5)) < tol & abs(planes(:,2) - planes(:,8)) < tol);
% horizontalInds = find(abs(planes(:,3) - planes(:,6)) < tol & abs(planes(:,3) - planes(:,9)) < tol);
% walls = union(verticalInds1,verticalInds2);
% walls = union(walls,horizontalInds);
% allInds = 1:1:length(faces);
% allInds(walls) =[];
% materialSurfs = allInds;
% materialZ(allInds) = 74.0;


% plotSet = walls;
plotSet = 1:1:length(planes);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

[X,Y,Z] = refineXYZ(X,Y,Z,6)
planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(10)
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','k')%'none')
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);

title({'Proto-MPEX Simulated GITR Geometry',  'for W Isotope Erosion Experiments'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('W 172','W 173','W 174','W 175')
% plotSet = allInds;
% surfs(allInds) = 1.0;
% 
% hold on
% X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
% Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
% Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
% figure(10)
% patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',.3,'EdgeColor','k')%'none')

nFaces = length(planes);
abcd = zeros(length(planes),4);
area = zeros(nFaces,1);
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
end
% X = [planes(:,1),planes(:,4),planes(:,7)];
% Y = [planes(:,2),planes(:,5),planes(:,8)];
% Z = [planes(:,3),planes(:,6),planes(:,9)];
% patch(transpose(X),transpose(Y),transpose(Z),'green','FaceAlpha',.3)
inDir = ones(1,length(planes));
figure(10)
hold on
for i=1:length(surfs)
    ind = i;
    normal = -abcd(ind,1:3)/plane_norm(ind);
    l_normal = 0.01;
    normal = l_normal*normal;
    centroid = [1/3*(planes(ind,1)+planes(ind,4)+planes(ind,7)), ...
                1/3*(planes(ind,2)+planes(ind,5)+planes(ind,8)), ...
                1/3*(planes(ind,3)+planes(ind,6)+planes(ind,9))];
        if(centroid(1) > 0.005 & normal(1)>0 & normal(2) == 0 & normal(3) == 0)
        inDir(ind) = -1;
        end    
        if(centroid(1) < -0.005 & normal(1)<0 & normal(2) == 0 & normal(3) == 0)
        inDir(ind) = -1;
        end 
        if(centroid(2) > 0.005 & normal(2)>0 & normal(1) == 0 & normal(3) == 0)
        inDir(ind) = -1;
        end    
        if(centroid(2) < -0.005 & normal(2)<0 & normal(1) == 0 & normal(3) == 0)
        inDir(ind) = -1;
        end   
        if(centroid(3) > 0.02 & normal(3)>0 & normal(2) == 0 & normal(1) == 0)
        inDir(ind) = -1;
        end    
        if(centroid(3) < 0.02 & normal(3)<0 & abs(normal(2))<1e-10)
        inDir(ind) = -1;
        end   
normal = inDir(ind)*normal;
    quiver3(centroid(1),centroid(2),centroid(3),normal(1),normal(2),normal(3))
end

fileID = fopen('gitrGeometryPointPlane3d.cfg','w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,'%5e',planes(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,1));
end
fprintf(fileID,' ] \n   y1 = [');
fprintf(fileID,'%5e',planes(1,2));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,2));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,'%5e',planes(1,3));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,3));
end
fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,'%5e',planes(1,4));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,4));
end
fprintf(fileID,' ] \n   y2 = [');
fprintf(fileID,'%5e',planes(1,5));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,5));
end
fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,'%5e',planes(1,6));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,6));
end
fprintf(fileID,' ] \n   x3 = [');
fprintf(fileID,'%5e',planes(1,7));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,7));
end
fprintf(fileID,' ] \n   y3 = [');
fprintf(fileID,'%5e',planes(1,8));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,8));
end
fprintf(fileID,' ] \n   z3 = [');
fprintf(fileID,'%5e',planes(1,9));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,9));
end
fprintf(fileID,' ] \n   a = [');
fprintf(fileID,'%5e',abcd(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,1));
end

fprintf(fileID,' ] \n   b = [');
fprintf(fileID,'%5e',abcd(1,2));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,2));
end

fprintf(fileID,' ] \n   c = [');
fprintf(fileID,'%5e',abcd(1,3));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,3));
end
fprintf(fileID,' ] \n   d = [');
fprintf(fileID,'%5e',abcd(1,4));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,4));
end
fprintf(fileID,' ] \n   plane_norm = [');
fprintf(fileID,'%5e',plane_norm(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',plane_norm(i));
end
% fprintf(fileID,' ] \n   ABxAC = [');
% fprintf(fileID,'%5e',ABxAC(1))
% for i=2:nFaces
% fprintf(fileID, ',')
% fprintf(fileID,'%5e',ABxAC(i))
% end
fprintf(fileID,' ] \n   BCxBA = [');
fprintf(fileID,'%5e',BCxBA(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',BCxBA(i));
end
fprintf(fileID,' ] \n   CAxCB = [');
fprintf(fileID,'%5e',CAxCB(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',CAxCB(i));
end
fprintf(fileID,' ] \n   area = [');
fprintf(fileID,'%5e',area(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',area(i));
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