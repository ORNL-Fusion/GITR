clear all
close all
load('D3DrzZ2D2Tiles10Res.mat')
figure(1)
plot(r,z)
hold on
scatter(r,z)
axis equal
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c);

surface = zeros(1,length(Z));
surface(find(Z>0))=1;

lines = GITR_LinesFromPoints(r,z, 'closed');
% Z = ones(1,length(r))*6.0;
GITR_writeGeomCFG_fromLines(lines,surface,'gitrD3DGeometry2DWrings.cfg',Z,0,0,0)