close all
clear all

file = 'positions.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
hitWall =ncread(file,'hitWall');
nP=length(x);
hit = find(hitWall);
notHit = find(hitWall==0);
He_Wspyl = [4.07e-3 3.26e-3  3.02e-3];
He_Wroughness = [1.0 1.1 1.2];
E = 0.5*12*1.66e-27*(vx.^2 + vy.^2 + vz.^2)/1.602e-19;
vvv = sign(vz).*(vx.^2 + vy.^2 + vz.^2);
Ehit = E(hit);
chargeHit = charge(hit);
r = sqrt(x.^2 + y.^2);
foundIndices = find(r<0.005);
figure(1)

scatter3(x,y,z)
hold on
%scatter3(x(hit),y(hit),z(hit))
title('Final W Particle Positions')
xlabel('X')
ylabel('Y')
zlabel('Z')

specFile = 'spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
% gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');
figure(2)
h=pcolor(gridR,gridZ,dens(:,:,5)')
h.EdgeColor = 'none';
figure(21)
semilogy(gridZ,sum(dens(:,:,5),1))
hold on
histogram(z(notHit))
maxDens = max(max(dens(:,:,5)));
 figure(21)
hold on
s = linspace(0,1.5);
plot(s,1e5*exp(-1*s))
close all
T=50;
m=12;
vTh = sqrt(T*1.602e-19/m/1.66e-27);
k = 1.38e-23*11604; 
B = m*1.66e-27/(2*T*k);
vgrid = linspace(-5*vTh,5*vTh);
fv1 = sqrt(B/pi)*exp(-B*vgrid.^2);
figure(1)


vtot = sign(vz).*sqrt(vz.*vz + vx.*vx + vy.*vy);
figure(2)
h1 = histogram(vz,[vgrid, vgrid(end) + (vgrid(2)-vgrid(1))]-0.5*(vgrid(2)-vgrid(1)));
h1binAve = mean([h1.BinEdges(1:end-1);h1.BinEdges(2:end)]);
h1Vals = h1.Values;
hold on
semilogy(vgrid,1.05*fv1*max(h1Vals)/max(fv1))
title({'Impurity Velocity (z-component) Distribution','In D plasma T=50 eV, n=1e19 m-3'})
xlabel('v_z [m/s]')
ylabel('Counts')
set(gca,'fontsize',16)
axis([vgrid(1) vgrid(end) 0 max(h1Vals)])
figure(1)
hold on
semilogy(sign(h1binAve).*(h1binAve./vTh).^2,h1Vals./sum(h1Vals))
semilogy(sign(vgrid).*(vgrid./vTh).^2,2000*fv1)
set(gca, 'YScale', 'log')
% axis([-(vgrid(1)/vTh)^2 (vgrid(1)/vTh)^2 min(h1Vals) 1.2*max(h1Vals)])
title({'Impurity Velocity (z-component) Distribution','In D plasma T=50 eV, n=1e19 m-3'})
xlabel('(v_z/v_{Th})^2')
ylabel('Counts')
set(gca,'fontsize',16)
file = 'history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
E = 0.5*12*1.66e-27*(vx.^2 + vy.^2 + vz.^2)/1.602e-19;
figure(19)
hold on
subX=5;
subY=5;
size_hist = size(vz);
time_snap = 1:(size_hist(1)-1)/subX/subY:size_hist(1);
% time_snap = 1:1:25;
dt = 1e-7;
for i=1:subY
    for j=1:subX
subplot(subX,subY,subX*(i-1) + j)
histogram(vz(time_snap(subX*(i-1) + j),:))
title(['t=', num2str(time_snap(subX*(i-1) + j)*dt),'s'])
axis([-1e5 1e5 0 10000])
    end

end
figure(20)
Eend = 0.5*12*1.66e-27*(vx(end,:).^2 + vy(end,:).^2 + vz(end,:).^2)/1.602e-19;
plot(mean(E'))
sizeArray = size(x);
nP = sizeArray(2);
% nT = 1e2;
% nP = 1e3;
% 
% x = reshape(x,[nP,nT]);
% y = reshape(y,[nP,nT]);
% z = reshape(z,[nP,nT]);
% charge = reshape(charge,[nP,nT]);
% figure(15)
% for i=1:nP
% plot3(x(i,:),y(i,:),z(i,:))
% hold on
% end

r = sqrt(x.^2 + y.^2);
figure(16)
for i=1:100
plot3(x(:,i),y(:,i),z(:,i))
hold on
end
ave_vz=mean(vz');
figure(101)
plot(ave_vz)
figure(102)
histogram(vz(end,:))
hold on
histogram(vx(end,:))
histogram(vy(end,:))
% n = reshape(dens(:,:,5),[length(gridR),length(gridZ)]);
% n = reshape(sum(n,2),[length(gridR),length(gridZ)]);
% figure(17)
% p3 = pcolor(gridR,gridZ,n')
% hold on
% plot(gridR,gridR*sind(30))

% figure(7)
%  [h1,c1] = hist3([x(hit) y(hit)],[30 30])
% % set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% % colorbar
% p1=pcolor(c1{1},c1{2},log10(h1'))
% colorbar
% title({'GITR W Target He Erosion Distribution','He/D Plasma Exposure 180205'})
% xlabel('x [m]') % x-axis label
% ylabel('y [m]') % y-axis label
% zlabel('z [m]') % y-axis label
% set(gca,'fontsize',16)
% % S = findobj('type','surf');
% % ZD = get(S,'zdata');
% % ZD(~ZD) = .1;
% % set(S,'zdata',ZD);
% % set(gca,'zscale', 'log');
% xs = x(hit);
% ys = y(hit);
% xg = c1{1};
% dx = xg(2)-xg(1);
% yg = c1{2};
% dy = yg(2)-yg(1);
% energies = zeros(length(xg),length(yg));
% charges = zeros(length(xg),length(yg));
% for i=1:length(xg)
% fx = find((xs >= xg(i)-0.5*dx) & (xs < xg(i)+0.5*dx));
%     for j=1:length(yg)
%     fy = find((ys >= yg(j)-0.5*dy) & (ys < yg(j)+0.5*dy));
%     fs = intersect(fx,fy);
% %     mean(Ehit(fs))
%     energies(i,j) = mean(Ehit(fs));
%     charges(i,j) = mean(chargeHit(fs));
%     end
% end
% 
% figure(8)
% p2=pcolor(c1{1},c1{2},energies')
% colorbar
% title({'GITR W Target He Erosion Distribution','He/D Plasma Exposure 180205'})
% xlabel('x [m]') % x-axis label
% ylabel('y [m]') % y-axis label
% zlabel('z [m]') % y-axis label
% set(gca,'fontsize',16)
% 
% figure(9)
% p2=pcolor(c1{1},c1{2},charges')
% colorbar
% title({'GITR W Target He Erosion Distribution','He/D Plasma Exposure 180205'})
% xlabel('x [m]') % x-axis label
% ylabel('y [m]') % y-axis label
% zlabel('z [m]') % y-axis label
% set(gca,'fontsize',16)
maxDens = max(max(dens(:,:,5)));
Abin = (gridR(2)-gridR(1))*(gridZ(2)-gridZ(1));
A0 = 0.1*0.04;
dt = 5e-8;
np = maxDens*A0*dt/nP/Abin
file = 'surface.nc';
grossDep = ncread(file,'grossDeposition');
grossEro = ncread(file,'grossErosion');
edist = ncread(file,'surfEDist');
spylDist = ncread(file,'surfSputtDist');
rfylDist = ncread(file,'surfReflDist');
sumWeightStrike = ncread(file,'sumWeightStrike');
aveSpyl = ncread(file,'aveSpyl');
spylCounts = ncread(file,'spylCounts');
sumParticlesStrike = ncread(file,'sumParticlesStrike');

grossDeposition = sum(grossDep)
grossErosion = sum(grossEro)
netEro = grossErosion - grossDeposition - nP
averageY = sum(aveSpyl)/sum(spylCounts)
sum(sumParticlesStrike)


eDtotal = reshape(sum(edist,3),[30,200]);
figure(18)
p4 = pcolor(linspace(0,100,200),linspace(0,90,30),eDtotal);
eaDist = sum(edist,3);
sDist = sum(spylDist,3);
rDist = sum(rfylDist,3);
% eaDist(eaDist > 2e3) = 0;
eaDist1 = reshape(eaDist,30,200);
sputtDist1 = reshape(sDist,30,200);
reflDist1 = reshape(rDist,30,200);
energy = linspace(0,99.5,200);
angle = linspace(0,87,30);
figure(20)
p1 = pcolor(energy,angle,sputtDist1)
p1.EdgeColor = 'none';
colorbar
title({'ITER Full Power Operation', 'Summed Target W Ion Energy-Angle Distribution'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label
figure(21)
p1 = pcolor(energy,angle,reflDist1)
p1.EdgeColor = 'none';
colorbar
title({'ITER Full Power Operation', 'Summed Target W Ion Energy-Angle Distribution'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label

if (exist('x1') == 0)
    fid = fopen('../input/gitrGeometryRev.cfg');
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:11
        tline = fgetl(fid);
        evalc(tline);
    end
length_line = length;
surfaces=surface;
clear length
clear surface
end

figure(19)
plot([x1 x1(1)],[z1 z1(1)])
x = x1;
y = z1;
z = zeros(size(x));
col = [grossDep; 0]';  % This is the color, vary with x in this case.

hold on

surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',4);
colorbar


file = 'history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
sizeArray = size(x);
nP = sizeArray(2);
% nT = 1e2;
% nP = 1e3;
% 
% x = reshape(x,[nP,nT]);
% y = reshape(y,[nP,nT]);
% z = reshape(z,[nP,nT]);
% charge = reshape(charge,[nP,nT]);
% figure(15)
% for i=1:nP
% plot3(x(i,:),y(i,:),z(i,:))
% hold on
% end

r = sqrt(x.^2 + y.^2);
figure(16)
for i=1:100
plot3(x(:,i),y(:,i),z(:,i))
hold on
end
sum_path = 0;
% for i=1:nP
%     whereIonized = find(charge(i,:));
%     if length(whereIonized) > 0
%         thisPart = sqrt((x(i,whereIonized(1))-x(i,1))^2 + (y(i,whereIonized(1))-y(i,1))^2 ...
%             + (z(i,whereIonized(1))-z(i,1))^2);
%         sum_path = sum_path + thisPart;
%     end
% end

% sum_path = sum_path/nP

specFile = 'spec.nc'
dens = ncread(specFile,'n');
gridR = ncread(specFile,'gridR');
gridY = ncread(specFile,'gridY');
gridZ = ncread(specFile,'gridZ');

n = reshape(dens(:,:,:,7),[length(gridR),length(gridY),length(gridZ)]);
n = reshape(sum(n,2),[length(gridR),length(gridZ)]);
figure(17)
p3 = pcolor(gridR,gridZ,n')
hold on
plot(gridR,gridR*sind(30))

file = 'surface.nc';
grossDep = ncread(file,'grossDeposition');
grossEro = ncread(file,'grossErosion');
edist = ncread(file,'surfEDist');
sumWeightStrike = ncread(file,'sumWeightStrike');
aveSpyl = ncread(file,'aveSpyl');
spylCounts = ncread(file,'spylCounts');
sumParticlesStrike = ncread(file,'sumParticlesStrike');
eDtotal = reshape(sum(edist,3),[30,200]);
figure(18)
p4 = pcolor(linspace(0,1000,200),linspace(0,90,30),eDtotal);

if (exist('x1') == 0)
    fid = fopen('../input/gitrGeometryRev.cfg');
    
    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:11
        tline = fgetl(fid);
        evalc(tline);
    end
length_line = length;
clear length
end

figure(19)
plot([x1 x1(1)],[z1 z1(1)])
x = x1;
y = z1;
z = zeros(size(x));
col = [grossDep; 0]';  % This is the color, vary with x in this case.

hold on

surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',4);
colorbar