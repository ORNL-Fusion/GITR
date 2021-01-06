
clearvars -except Data
%clear all
close all
clc
load('Mask.mat')

OxyPercent=0.01;

%erosion_rate=2.6988e18; %Low density Non-Mag Case
erosion_rate=5.8723e21; %Low density Mag Case
erosion_rate=6e18;
erosion_rate = 3.4044e+18;
%erosion_rate=5.18e20; %High density Mag Case
%erosion_rate=0; High density Non-Mag Case

%erosion_rate=3.9036e18; %Low density High Te Non-Mag Case
%erosion_rate=8.434018812449574e+21; %Low density High Te Mag Case
%erosion_rate=8.434018812449574e+21; %Low density High Te 1mill particles Case
%erosion_rate=3.091855020025960e+22; %Low density High Te Finer Mesh Case
%erosion_rate=7.421398052258447e+20; %High density High Te Mag Case

%erosion_rate=2.661066105959504e+21; %Low density High Te Non-Mag Oxygen Case
%erosion_rate=4.725712484384497e+22; %Low density High Te Mag Oxygen Case
%erosion_rate=1.916989293644594e+20; %High density High Te Non-Mag Oxygen Case
%erosion_rate=1.263748583934766e+22; %High density High Te Mag Oxygen Case


%erosion_rate=2.8727e20; %100V LowDensity HighTe

%erosion_rate=2.001729188974246e+20; %100V LowDensity
%erosion_rate=1.065297151921256e+21; %100V HighDensity
%erosion_rate=8.717812559255756e+20; %500V LowDensity
%erosion_rate=4.596827001220502e+21; %500V HighDensity
%erosion_rate=1.036316162578123e+21; %1kV LowDensity
%erosion_rate=5.461506751394209e+21; %1kV HighDensity

nP=100000;
erosionPP=erosion_rate/nP;

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
    
plot_tracks = 0;
file = strcat(pwd,'/output/positions.nc');
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
file = strcat(pwd,'/output/history.nc');
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
    fid = fopen(strcat(pwd,'/gitrGeometryPointPlane3d.cfg'));
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end

%%
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
set(gca,'fontsize',13)
set(gcf,'color','w')
% set(cb1,'fontsize',14);
axis equal
%%

file = strcat(pwd,'/output/surface.nc');
grossDep0 = ncread(file,'grossDeposition');
grossEro0 = ncread(file,'grossErosion');

subset = find(Zsurface);
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
% THETA = atan2d(Y,X);
% ss=sum(sign(THETA),2);
% swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
% THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
THETA = 0*Z;
for ii=1:length(THETA)
    %     if ii==2457
%     if ii==4859
%         ii
%     end
    THETA(ii,:) = atan2d(Y(ii,:),X(ii,:));
    diff21 = abs(THETA(ii,2) - THETA(ii,1));
    diff31 = abs(THETA(ii,3) - THETA(ii,1));
    diff32 = abs(THETA(ii,3) - THETA(ii,2));
    
    maxdiff = max([diff21,diff31,diff32]);
    if maxdiff > 90
        if THETA(ii,1)>0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)<0
                THETA(ii,2)=THETA(ii,2)+360;
            end
            if THETA(ii,3)<0
                THETA(ii,3)=THETA(ii,3)+360;
            end
        end
        if THETA(ii,1)<0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)>0
                THETA(ii,2)=THETA(ii,2)-360;
            end
            if THETA(ii,3)>0
                THETA(ii,3)=THETA(ii,3)-360;
            end
        end
    end
end
THETA=THETA+180;

%%
ElementArea=area(surface)';
pointsize=10;

GrossErosion=(erosionPP.*grossEro0./ElementArea);
%GrossErosion=GrossErosion*OxyPercent+((1-OxyPercent)*Data.GrossErosion);

figure(2)
patch(transpose(Z),transpose(THETA),0*transpose(Z),GrossErosion,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
scatter(Mask(:,1),Mask(:,2),pointsize,(max(GrossErosion)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
% patch(transpose(X),transpose(Y),transpose(Z),grossEro0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
c=colorbar;
ylabel(c, 'Particles/m^2/s','FontSize',13)
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
title('Gross Erosion','FontSize',13);
hold off
% set(gca,'ColorScale','log')
% caxis([1e22 1e23])

%view(90,90)
%%

GrossDeposition=erosionPP.*grossDep0./ElementArea;
%GrossDeposition=GrossDeposition*OxyPercent+((1-OxyPercent)*Data.GrossErosion);

figure(3)
patch(transpose(Z),transpose(THETA),0*transpose(Z),GrossDeposition,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
scatter(Mask(:,1),Mask(:,2),pointsize,(max(GrossDeposition)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

% patch(transpose(X),transpose(Y),transpose(Z),grossDep0,'FaceAlpha',.3,'EdgeAlpha', 0.3)
c=colorbar;
ylabel(c, 'Particles/m^2/s','FontSize',13)
title('Net Erosion','FontSize',13);
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
caxis([0 inf])
set(gcf,'color','w')
title('Gross Deposition','FontSize',13);
hold off

% THETA = atan2d(Y,X);
% ss=sum(sign(THETA),2);
% swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
% THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
% 
% for ii=1:length(THETA)
%     if THETA(ii,1)<0
% THETA(ii,1)=THETA(ii,1)+360;
%     end
%     if THETA(ii,2)<0
% THETA(ii,2)=THETA(ii,2)+360;
%     end
%     if THETA(ii,3)<0
% THETA(ii,3)=THETA(ii,3)+360;
%     end
% end

NetErosion=erosionPP.*(grossEro0-grossDep0)./ElementArea;
%NetErosion=NetErosion*OxyPercent+((1-OxyPercent)*Data.NetErosion);

figure(4)
patch(transpose(Z),transpose(THETA),0*transpose(Z),NetErosion,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
scatter(Mask(:,1),Mask(:,2),pointsize,(max(NetErosion)*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel(c, 'Particles/m^2/s','FontSize',13)
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
caxis([0 inf])
set(gcf,'color','w')
title('Net Erosion of Al','FontSize',13);
hold off

netEro = (grossEro0-grossDep0);
netEro(find(netEro >= 0)) = nan;

% figure(5)
% patch(transpose(Z),transpose(THETA),0*transpose(Z),erosionPP.*netEro./ElementArea,'FaceAlpha',1,'EdgeAlpha', 0.3)
% hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netEro./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);
% 
% c=colorbar;
% ylabel(c, 'Particles/m^2/s','FontSize',13);
% ylabel('Angle [Deg]','FontSize',13);
% xlabel('Z [m]','FontSize',13);
% set(gcf,'color','w')
% caxis([0 inf])
% title({'Net Erosion of Al','Deposition Regions Only'},'FontSize',13);
% hold off

netDep=abs(netEro);

NetDep=erosionPP.*netDep./ElementArea;
%NetDep=NetDep*OxyPercent+((1-OxyPercent)*Data.NetDeposition;

for ii=1:length(NetDep)
   if isnan(NetDep(ii))==1
       NetDep(ii)=0;
   end
end
Z_cyl = Z;
figure(6)
patch(transpose(Z),transpose(THETA),0*transpose(Z),NetDep,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel(c, 'Particles/m^2/s','FontSize',13);
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
title({'Net Deposition of Al'},'FontSize',13);
hold off

%% End cap physics
nExited = sum(weight0(find(abs(z0)>=0.2)));

ExitedParticles=nExited*erosionPP;
PercentExited=nExited/nP*100;

%%% histogram of particles which hit end cap
nBins = 20;
edges = linspace(-0.07,0.07,nBins+1);
centers = linspace(-0.07,0.07,nBins);
hist_bins = zeros(1,nBins);

for i=1:nBins

    hist_bins(i) = sum(weight0(find(z0>=0.1999 & x0 > edges(i) & x0 <= edges(i+1))));

end

figure(321)
plot(centers,hist_bins)
xlabel('X [cm]','FontSize',13);
ylabel('Particle Count','FontSize',13);
set(gcf,'color','w')
title({'Paritlces leaving the helicon'},'FontSize',13);

%%

%For total 
% hitCap1 = find(z0<=-0.199999);
% hitCap2 = find(z0>=0.1999);

%For nuetrals 
% hitCap1 = find(z0<=-0.199999 & charge0==0);
% hitCap2 = find(z0>=0.1999 & charge0==0);

%For ionized
hitCap1 = find(z0<=-0.1999 & charge0>=0);
hitCap2 = find(z0>=0.1999 & charge0>=0);

figure(322)
s1=scatter3(x0(hitCap1),y0(hitCap1),z0(hitCap1));
hold on
s2=scatter3(x0(hitCap2),y0(hitCap2),z0(hitCap2));
xlabel('X [cm]','FontSize',13);
ylabel('Y [cm]','FontSize',13);
zlabel('Z [cm]','FontSize',13);
set(gcf,'color','w')
title({'Paritlces leaving the helicon'},'FontSize',13);
hold off

% end cap flux 

%Upstream based on geometry
cap1 = find(Zsurface ==0 & z1 < 0);

cap1_bins = zeros(1,length(cap1));
for i=1:length(cap1)
  in = inpolygon(x0(hitCap1),y0(hitCap1),[x1(cap1(i)),x2(cap1(i)),x3(cap1(i))],[y1(cap1(i)),y2(cap1(i)),y3(cap1(i))]);
  cap1_bins(i) = sum(weight0(hitCap1(in)));
end
subset = cap1;%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(333)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),erosionPP*cap1_bins./area(cap1),'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Upstream End Cap Flux ','FontSize',13);
xlabel('X [m]','FontSize',13);
ylabel('Y [m]','FontSize',13);
zlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gca,'fontsize',13)
% set(cb1,'fontsize',14);
axis equal
c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',13)

%Downstream based on geometry
cap2 = find(Zsurface ==0 & z1 > 0);
cap2_bins = zeros(1,length(cap2));
for i=1:length(cap2)
  in = inpolygon(x0(hitCap2),y0(hitCap2),[x1(cap2(i)),x2(cap2(i)),x3(cap2(i))],[y1(cap2(i)),y2(cap2(i)),y3(cap2(i))]);
  cap2_bins(i) = sum(weight0(hitCap2(in)));
end
subset = cap2;%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(334)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];

DownstreamFlux=erosionPP*cap2_bins./area(cap2);
%DownstreamFlux=DownstreamFlux*OxyPercent+((1-OxyPercent)*Data.DownstreamFlux);

%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),DownstreamFlux,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Downstream End Cap Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% axis([xMin xMax yMin yMax zMin zMax])
% cb1 = colorbar
set(gcf,'color','w')
set(gca,'fontsize',13)
% set(cb1,'fontsize',14);
axis equal
c=colorbar;
ylabel(c, 'Flux [Particles/m^2/s]','FontSize',13)

% figure
% histogram(charge0(hitCap1))
% figure
% histogram(charge0(hitCap2))



%%
% figure
% s3=scatter(x0(hitCap2),y0(hitCap2));
% s3.MarkerFaceAlpha =weight0(hitCap1);

%%

Data.GrossErosion=(erosionPP.*grossEro0./ElementArea);
Data.GrossDeposition=erosionPP.*grossDep0./ElementArea;
Data.NetErosion=erosionPP.*(grossEro0-grossDep0)./ElementArea;
Data.NetDeposition=erosionPP.*netDep./ElementArea;
Data.DownstreamFlux=erosionPP*cap2_bins./area(cap2);

%%
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

aveSpyl = ncread(file,'aveSpyl');
spylCounts = cast(ncread(file,'spylCounts'),'single');
sumWeightStrike = ncread(file,'sumWeightStrike');
surfSputtDist = ncread(file,'surfSputtDist');

eDist = ncread(file,'surfEDist');
e_grid = 10:10:1000;
s1=reshape(sum(eDist,1),100,14744);
s2 = s1.*e_grid'./sum(s1);
meanE = sum(s2);

sp1=reshape(sum(surfSputtDist,1),100,14744);
sp2 = sp1.*e_grid'./sum(sp1);
mean_spE = sum(sp2);

yield =aveSpyl./spylCounts;

figure(201)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),meanE,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
title({'meanE'},'FontSize',13);
hold off

figure(202)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),yield,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
title({'yield'},'FontSize',13);
hold off

figure(203)
patch(transpose(Z_cyl),transpose(THETA),0*transpose(Z_cyl),mean_spE,'FaceAlpha',1,'EdgeAlpha', 0.3)
hold on
% scatter(Mask(:,1),Mask(:,2),pointsize,(max((erosionPP.*netDep./ElementArea))*Mask(:,3)),'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);

c=colorbar;
ylabel('Angle [Deg]','FontSize',13);
xlabel('Z [m]','FontSize',13);
set(gcf,'color','w')
title({'mean sputtered energy'},'FontSize',13);
hold off
