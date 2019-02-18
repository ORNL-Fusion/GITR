function GITR_writeGeomCFG_fromLines(lines,surfaces,filename,Z,y1,y2,periodic)
nPoints = length(lines(:,1));

fileID = fopen(filename,'w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,'%f',lines(1,1));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,1));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,'%f',lines(1,2));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,2));
end

fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,'%f',lines(1,3));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,3));
end

fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,'%f',lines(1,4));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,4));
end

fprintf(fileID,' ] \n   slope = [');
fprintf(fileID,'%f',lines(1,5));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,5));
end

fprintf(fileID,' ] \n   intercept = [');
fprintf(fileID,'%f',lines(1,6));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,6));
end

fprintf(fileID,' ] \n   length = [');
fprintf(fileID,'%f',lines(1,7));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',lines(i,7));
end

fprintf(fileID,' ] \n   Z = [');
fprintf(fileID,'%f',Z(1));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%f',Z(i));
end
fprintf(fileID, ',');
fprintf(fileID,'%f',0.0);
fprintf(fileID, ',');
fprintf(fileID,'%f',0.0);
fprintf(fileID,' ] \n   surface = [');
fprintf(fileID,'%i',surfaces(1));
for i=2:nPoints
fprintf(fileID, ',');
fprintf(fileID,'%i',surfaces(i));
end

fprintf(fileID,' ] \n   y1 = ');
fprintf(fileID,'%f',y1);

fprintf(fileID,' ; \n   y2 = ');
fprintf(fileID,'%f',y2);

fprintf(fileID,' ; \n   periodic = ');
fprintf(fileID,'%i',periodic);

fprintf(fileID,' ; \n}');