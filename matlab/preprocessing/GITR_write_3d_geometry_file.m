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

figure(10)
hold on
Wsurfs = find(materialZ);
inDir = ones(1,length(planes));

for i=1:length(Wsurfs)
    ind = Wsurfs(i);
    normal = -abcd(ind,1:3)/plane_norm(ind);
    if(normal(3) <0)
        inDir(ind) = -1;
    end
    quiver3(planes(ind,1),planes(ind,2),planes(ind,3),normal(1),normal(2),inDir(ind)*normal(3))
end
fileID = fopen(geometryFilename,'w');
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
fprintf(fileID,'theta1 = 0.0\n }');
fclose(fileID)