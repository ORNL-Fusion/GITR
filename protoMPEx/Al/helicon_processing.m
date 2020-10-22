clear all
close all
% M = csvread('iterGeom.csv');
% r = M(:,1);
% z = M(:,2);
% rGeom = r;
% zGeom = z;
total_redeposition_fraction = [];
local_redep_fraction = [];
prompt_redep_fraction = [];
net_erosion_fraction = [];
angles = 0:10:90;
% for angle=0:10:90
    
plot_tracks = 1;
file = 'output/positions.nc';
hitWall = ncread(file,'hitWall');
nHit = length(find(hitWall));
hasHit = find(hitWall);
notHit = find(hitWall==0);
x0 = ncread(file,'x');
y0 = ncread(file,'y');
z0 = ncread(file,'z');
vx0 = ncread(file,'vx');
vy0 = ncread(file,'vy');
vz0 = ncread(file,'vz');
distTraveled = ncread(file,'distTraveled');
charge0 = ncread(file,'charge');
weight0 = ncread(file,'weight');
vtot = sqrt(vx0.^2 +vy0.^2 + vz0.^2);
E = 0.5*27*1.66e-27*vtot.^2/1.602e-19;
figure(11)
histogram(E)

figure(10)
scatter3(x0(hasHit),y0(hasHit),z0(hasHit))
hold on
scatter3(x0(notHit),y0(notHit),z0(notHit))
local_length = 8e-4;
prompt_length = 31e-3;
redep = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1));
redep1 = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & (charge0 == 1));
redep2 = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & (charge0 == 2));
redep_local = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & distTraveled < local_length);
redep_prompt = find((z0 > -0.001) & (z0 < 0.00) & (hitWall == 1) & distTraveled > local_length & distTraveled < prompt_length);

if plot_tracks
file = 'output/history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
sizeArray = size(x);
nP = sizeArray(2);
hit = find(weight(end,:) < 1);

r = sqrt(x.^2 + y.^2);
 figure(1)
 hold on
%  scatter3(x0(redep1),y0(redep1),z0(redep1))
% % i=581
% % plot(r(:,i),z(:,i))
% %41242
for i=1:1:length(x0)
 
% plot(r(:,i),z(:,i))
plot3(x(:,i),y(:,i),z(:,i),'k')
end
end
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% % plot(rGeom,zGeom)
% % hold on
% % scatter(rGeom,zGeom)
% axis equal
hold on
if (exist('x1') == 0)
    fid = fopen('input/gitrGeometryHelicon.cfg');
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end

surface = find(Zsurface);
nSurfaces = length(a);
%Section for finding a subset of planes to plot
r = sqrt(x1.^2 + y1.^2);
subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface)
title('Sample Particle Tracks for Helicon Erosion Source')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gca,'fontsize',16)
% set(cb1,'fontsize',14);
axis equal

file = 'output/surface.nc';
grossDep0 = ncread(file,'grossDeposition');
grossEro0 = ncread(file,'grossErosion');

subset = find(Zsurface);
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
THETA = atan2(Y,X);
ss=sum(sign(THETA),2);
swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );

figure(2)
patch(transpose(THETA),transpose(Z),0*transpose(Z),grossEro0,'FaceAlpha',1,'EdgeAlpha', 0.3)

% patch(transpose(X),transpose(Y),transpose(Z),grossEro0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
colorbar
title('Net Erosion')
xlabel('theta [rad]')
ylabel('z [m]')
title('Gross Erosion')

figure(3)
patch(transpose(THETA),transpose(Z),0*transpose(Z),grossDep0,'FaceAlpha',1,'EdgeAlpha', 0.3)

% patch(transpose(X),transpose(Y),transpose(Z),grossDep0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
colorbar
title('Net Erosion')
xlabel('theta [rad]')
ylabel('z [m]')

title('Gross Deposition')

THETA = atan2(Y,X);
ss=sum(sign(THETA),2);
swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
figure(4)
patch(transpose(THETA),transpose(Z),0*transpose(Z),(grossEro0-grossDep0),'FaceAlpha',1,'EdgeAlpha', 0.3)
colorbar
title('Net Erosion')
xlabel('theta [rad]')
ylabel('z [m]')
% zlabel('z [m]')
title('Net Erosion of Al')

netEro = (grossEro0-grossDep0);
netEro(find(netEro >= 0)) = nan;
figure(5)
patch(transpose(THETA),transpose(Z),0*transpose(Z),netEro,'FaceAlpha',1,'EdgeAlpha', 0.3)
colorbar
title('Net Erosion')
xlabel('theta [rad]')
ylabel('z [m]')
% zlabel('z [m]')
title({'Net Erosion of Al','Deposition Regions Only'})

% figure(4)
% histogram(charge(end,:))

% figure(123)
% histogram(distTraveled(redep),0:0.0002:0.5)
% % hold on
% % histogram(distTraveled(redep2),0:0.0002:0.5)
% total_redeposition_fraction = [total_redeposition_fraction; length(redep)/length(x0)];
% local_redep_fraction = [local_redep_fraction; length(redep_local)/length(x0)];
% prompt_redep_fraction = [prompt_redep_fraction; length(redep_prompt)/length(x0) ];
% net_erosion_fraction = [net_erosion_fraction; (length(x0) - length(redep))/length(x0) ];
% % end
% close all
% figure(1)
% plot(angles,total_redeposition_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Total Re-deposited Fraction of Eroded W')
% axis([0 90 0.7 0.9])
% set(gca,'fontsize',16)
% 
% figure(2)
% plot(angles,local_redep_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Locally Re-deposited Fraction of Eroded W')
% axis([0 90 0 0.04])
% set(gca,'fontsize',16)
% 
% figure(3)
% plot(angles,prompt_redep_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Promptly Re-deposited Fraction of Eroded W')
% axis([0 90 0.25 0.7])
% set(gca,'fontsize',16)
% 
% figure(4)
% plot(angles,net_erosion_fraction,'-o')
% xlabel('Angle of Incidence [degrees]')
% ylabel('Fraction')
% title('Net Eroded Fraction of Eroded W')
% axis([0 90 0.1 0.25])
% set(gca,'fontsize',16)
