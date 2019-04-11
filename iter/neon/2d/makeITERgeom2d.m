clear all
close all
M = csvread('iterGeom.csv');
r = M(:,1);
z = M(:,2);
Z = zeros(1,length(r));
boundaryInd = [1:6,92:97,115:116];
Wind = [7:24,74:91,98:114];
BeInd= [25:73];
Z(boundaryInd) = 0;
Z(Wind) = 74;
Z(BeInd) = 4;

figure(1)
hold on
for i=1:length(boundaryInd)
    ind=boundaryInd(i);
    indPlus1 = ind+1;
    if indPlus1 > length(r)
        indPlus1=1;
    end
    p1=plot([r(ind),r(indPlus1)],[z(ind),z(indPlus1)],'k','lineWidth',2,'DisplayName','b')
end
for i=1:length(Wind)
    ind=Wind(i);
    indPlus1 = ind+1;
    if indPlus1 > length(r)
        indPlus1=1;
    end
    p2=plot([r(ind),r(indPlus1)],[z(ind),z(indPlus1)],'g','lineWidth',2,'DisplayName','W')
end
for i=1:length(BeInd)
    ind=BeInd(i);
    indPlus1 = ind+1;
    if indPlus1 > length(r)
        indPlus1=1;
    end
    p3=plot([r(ind),r(indPlus1)],[z(ind),z(indPlus1)],'b','lineWidth',2,'DisplayName','Be')
end
legend([p1;p2;p3])
axis equal
scatter(r,z)
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c);

surface = ones(1,length(Z));
% surface(find(Z>0))=1;

lines = GITR_LinesFromPoints(r,z, 'closed');
% Z = ones(1,length(r))*6.0;
GITR_writeGeomCFG_fromLines(lines,surface,'iterGeom2D.cfg',Z,0,0,0)