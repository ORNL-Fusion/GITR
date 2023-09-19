%cube geom
close all
clear all

% This is a manual specification of the triangles that make up a
% 3D surface geometry. CAD files and other sources can be used.

% The coordinates are specified in cartesion coordinates
% To give a cuboid of x [1,2], y [-1,1], and z [0,1]

% In the gitrGeometry.cfg file, first the point coordinates
% are specified. The three points that make up triangle i are
% A = [x1(i) y1(i) z1(i)], B = [x2(i) y2(i) z2(i)], and C = [x3(i) y3(i) z3(i)]

x1 = [1 1 1 1 2 2 1 1 1 1 1 1];
y1 = [-1 -1 1 1 -1 -1 -1 -1 -1 -1 -1 -1];
z1 = [0 0 0 0 0 0 0 0 1 1 0 0];

x2 = [2 1 2 1 2 2 1 1 2 1 2 1];
y2 = [-1 -1 1 1 1 -1 1 -1 -1 1 -1 1];
z2 = [0 1 0 1 0 1 0 1 1 1 0 0];

x3 = [2 2 2 2 2 2 1 1 2 2 2 2];
y3 = [-1 -1 1 1 1 1 1 1 1 1 1 1];
z3 = [1 1 1 1 1 1 1 1 1 1 0 0];

% Assembling the coordinates in this form is useful for matlab patch plots
X = [x1;x2;x3];
Y = [y1;y2;y3];
Z = [z1;z2;z3];

figure(11)
patch(X,Y,Z,'blue','FaceAlpha',0.5)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title("GITR Surface Mesh")
set(gca,'fontsize',15)

% This is a useful way of visualizing part of the gitr
% geometry, specifying individual or sets of indices to view
figure(12)
pind = 1; % Indices here. The find function is also useful here. E.g., find (z1 > 0.9)
patch(X(:,pind),Y(:,pind),Z(:,pind),'blue','FaceAlpha',0.5)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title("GITR Surface Mesh Subdomain")
set(gca,'fontsize',15)

% For many of the matlab scripts I will post, this is a useful format for
% the coordinates to be put in before additional processing.
planes = [x1' y1' z1' x2' y2' z2' x3' y3' z3'];
n_elements = length(x1);

% Here we define and initialize additional fields which are part of the
% GITR geometry
materialZ = ones(1,n_elements)*23; % Atomic number of the material surface


abcd = zeros(length(planes),4); % plane equation coefficients ax + by + cz = d
area = zeros(n_elements,1); % area in m^3
plane_norm = zeros(n_elements,1); % normalizing factor for the vector normal given by (a,b,c)
% BCxBA = zeros(n_elements,1); % pre-computed cross-products
% CAxCB = zeros(n_elements,1); % pre-computed cross-products


surfaces = ones(n_elements,1); % 1 or 0, collect data on this surface or not. This is for data management in large geometries

inDir = ones(n_elements,1); % 1 or -1 indicating if the default normal vector points
% into the volume or not. Normal vector = inDir*(a,b,c)/plane_norm.
% For simulations using the surface model, this requires some manual or
% programatic specification based on user knowledge
ind_flip = [1 4 5 8 9 12]; % Flip these normals to point into the geometry volume
inDir(ind_flip) = -1;

centroid = zeros(n_elements,3); % Only used for processing and plotting (quiver plots)

figure(11)
hold on
% Loop to calculate the needed fields
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
%     tmp = min(norm1,norm2);
%     shortestEdge(i) = min(tmp,norm3);
%     len = 0.01;
%     if (norm1 > len & norm2 > len & norm3> len)
%         allLong(i) = 1;
%     end

    s = (norm1+norm2+norm3)/2;
    area(i) = sqrt(s*(s-norm1)*(s-norm2)*(s-norm3));
    normalVec = cross(AB,AC);

    d = -(dot(normalVec,A));

    abcd(i,:) = [normalVec,d];
    plane_norm(i) = norm(normalVec);

%     BCxBA(i) = sign(dot(cross(BC,BA),normalVec));
%     CAxCB(i) = sign(dot(cross(CA,CB),normalVec));
    centroid(i,:) = [1/3*(planes(i,1)+planes(i,4)+planes(i,7)), ...
        1/3*(planes(i,2)+planes(i,5)+planes(i,8)), ...
        1/3*(planes(i,3)+planes(i,6)+planes(i,9))];
    quiver3(centroid(i,1),centroid(i,2),centroid(i,3),inDir(i)*abcd(i,1)/plane_norm(i),inDir(i)*abcd(i,2)/plane_norm(i),inDir(i)*abcd(i,3)/plane_norm(i),0.5)
end

format_ = '%5e';
fileID = fopen('gitrGeometry.cfg','w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,format_,planes(1,1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,1));
end
fprintf(fileID,' ] \n   y1 = [');
fprintf(fileID,format_,planes(1,2));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,2));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,format_,planes(1,3));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,3));
end
fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,format_,planes(1,4));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,4));
end
fprintf(fileID,' ] \n   y2 = [');
fprintf(fileID,format_,planes(1,5));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,5));
end
fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,format_,planes(1,6));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,6));
end
fprintf(fileID,' ] \n   x3 = [');
fprintf(fileID,format_,planes(1,7));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,7));
end
fprintf(fileID,' ] \n   y3 = [');
fprintf(fileID,format_,planes(1,8));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,8));
end
fprintf(fileID,' ] \n   z3 = [');
fprintf(fileID,format_,planes(1,9));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,planes(i,9));
end
fprintf(fileID,' ] \n   a = [');
fprintf(fileID,format_,abcd(1,1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,abcd(i,1));
end

fprintf(fileID,' ] \n   b = [');
fprintf(fileID,format_,abcd(1,2));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,abcd(i,2));
end

fprintf(fileID,' ] \n   c = [');
fprintf(fileID,format_,abcd(1,3));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,abcd(i,3));
end
fprintf(fileID,' ] \n   d = [');
fprintf(fileID,format_,abcd(1,4));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,abcd(i,4));
end
fprintf(fileID,' ] \n   plane_norm = [');
fprintf(fileID,format_,plane_norm(1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,plane_norm(i));
end
% fprintf(fileID,' ] \n   ABxAC = [');
% fprintf(fileID,'%5e',ABxAC(1))
% for i=2:n_elements
% fprintf(fileID, ',')
% fprintf(fileID,'%5e',ABxAC(i))
% end
% fprintf(fileID,' ] \n   BCxBA = [');
% fprintf(fileID,'%5e',BCxBA(1));
% for i=2:n_elements
%     fprintf(fileID, ',');
%     fprintf(fileID,'%5e',BCxBA(i));
% end
% fprintf(fileID,' ] \n   CAxCB = [');
% fprintf(fileID,'%5e',CAxCB(1));
% for i=2:n_elements
%     fprintf(fileID, ',');
%     fprintf(fileID,'%5e',CAxCB(i));
% end
fprintf(fileID,' ] \n   area = [');
fprintf(fileID,format_,area(1,1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,format_,area(i));
end
fprintf(fileID,' ] \n   Z = [');
fprintf(fileID,'%f',materialZ(1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,'%f',materialZ(i));
end
fprintf(fileID,' ] \n   surface = [');
fprintf(fileID,'%i',surfaces(1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,'%i',surfaces(i));
end
fprintf(fileID,' ] \n   inDir = [');
fprintf(fileID,'%i',inDir(1));
for i=2:n_elements
    fprintf(fileID, ',');
    fprintf(fileID,'%i',inDir(i));
end
fprintf(fileID,' ] \n}');
fprintf(fileID,'theta0=0.0;\n');
fprintf(fileID,'theta1=0.0;\n');
fprintf(fileID,'periodic=0;\n');
fclose(fileID);