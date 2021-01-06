clear all
close all
%Check Ftridyn file


% plasma_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_Plasma.txt';
% plasma_file = 'assets/Beers_Helicon_3D_100kW_HighDensity_Plasma.txt';
 plasma_file = 'assets/Beers_Helicon_3D_100kW_LowDensityHighTe_Plasma.txt';
% plasma_file = 'assets/Beers_Helicon_3D_100kW_HighDensityHighTe_Plasma.txt';

% surface_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle9.xlsx';
 surface_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle7_MagCase.xlsx';
% surface_file = 'assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle6.xlsx';
% surface_file = 'assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle11_MagCase.xlsx';
% surface_file ='assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle_100V.xlsx';
%surface_file ='assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle_100V.xlsx';
%surface_file ='assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle_500V.xlsx';
%surface_file ='assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle_500V.xlsx';
%surface_file ='assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle_1kV.xlsx';
%surface_file ='assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle_1kV.xlsx';

radius = 0.06256;
model = createpde;
importGeometry(model,'assets/helicon_6p256cm.stl');% Import STL file
figure(2)
pdegplot(model,'FaceLabels','on') %Plot stl 

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.005);% Options Hmax and Hmin can be set, linear order can also be used
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

% [X,Y,Z] = refineXYZ(X,Y,Z,6)
planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(10)
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','k')%'none')
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);

title({'Proto-MPEX Helicon Simulated GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('Al')
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
centroid = zeros(nFaces,3);
plane_norm = zeros(nFaces,1);
BCxBA = zeros(nFaces,1);
CAxCB = zeros(nFaces,1);
density = zeros(nFaces,1);

helicon_data

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
    density(i) = interpn(y,z,dens,rr,centroid(3));        
end
figure(11)
patch(transpose(X),transpose(Y),transpose(Z),density,'FaceAlpha',.3,'EdgeColor','k')%'none')
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);

title({'Proto-MPEX Simulated GITR Geometry',  'for W Isotope Erosion Experiments'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
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
%     centroid = [1/3*(planes(ind,1)+planes(ind,4)+planes(ind,7)), ...
%                 1/3*(planes(ind,2)+planes(ind,5)+planes(ind,8)), ...
%                 1/3*(planes(ind,3)+planes(ind,6)+planes(ind,9))];
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
figure(110)
patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',.3,'EdgeColor','k')%'none')


helicon_data

% [xx yy zz emag] = helicon_data_surface();
[xx yy zz emag ne_surf] = helicon_read_xls(surface_file);

for i=1:length(centroid)
    distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);
    
    [M I] = min(distance);
    e_value(i) = emag(I);
    dens_value(i) = ne_surf(I);
    if abs(centroid(i,3))> 0.15
        e_value(i) = 0;
        dens_value(i) = 0;
    end
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
fprintf(fileID,' ] \n   potential = [');
fprintf(fileID,'%5e',e_value(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%5e',e_value(i));
end
            fprintf(fileID,' ] \n');
            fprintf(fileID,'periodic = 0;\n');
            fprintf(fileID,'theta0 = 0.0;\n');
            fprintf(fileID,'theta1 = 0.0\n');
            fprintf(fileID,'periodic_bc_x0 = 0.0;\n');
            fprintf(fileID,'periodic_bc_x1 = 0.0;\n');
            fprintf(fileID,'periodic_bc_x = 0;}\n');
            fclose(fileID)
            


e_value = e_value; %Conversion of E-field to potential
plotSet = 1:1:length(planes);
plotSet = find(surfs>0);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
THETA = atan2(Y,X);
ss=sum(sign(THETA),2);
swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
figure(101)
% patch(transpose(X),transpose(Y),transpose(Z),e_value,'FaceAlpha',.3,'EdgeColor','k')
patch(transpose(THETA),transpose(Z),0*transpose(Z),e_value(plotSet),'FaceAlpha',1,'EdgeColor','k')
colorbar
title({'Surface Potential [V]','Low Flux Case'})
xlabel('Theta [radian]') % x-axis label
ylabel('z [m]') % y-axis label
set(gca,'fontsize',16)


file = 'assets/ftridynBackground.nc';
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, nE] = netcdf.inqDim(ncid,0);
[dimname, nA] = netcdf.inqDim(ncid,1);
if strcmp(file,'assets/ftridynBackground.nc')
[dimname, nS] = netcdf.inqDim(ncid,2);
else
    nS = 1;
end
energy = ncread(file,'E');
angle = ncread(file,'A');
spyld = ncread(file,'spyld');
rfyld = ncread(file,'rfyld');
cosxDist = ncread(file,'cosXDist');
cosxDistRef = ncread(file,'cosXDistRef');
cosyDist = ncread(file,'cosYDist');
% coszDist = ncread(file,'cosZDist');
eDist = ncread(file,'energyDist');
eDistRef = ncread(file,'energyDistRef');
eDistEgrid = ncread(file,'eDistEgrid');
eDistEgridRef = ncread(file,'eDistEgridRef');
phiGrid = ncread(file,'phiGrid');
thetaGrid = ncread(file,'thetaGrid');
thisEdistRef = reshape(eDistRef(:,1,:),length(eDistEgridRef),[]);
% figure(100)
% plot(eDistEgridRef,thisEdistRef)

spyld=reshape(spyld(:,:,1),length(angle),length(energy));

figure(113)
h = pcolor(energy,angle,spyld)
h.EdgeColor = 'none';
colorbar
set(gca,'ColorScale','log')
% set(gca, 'YDir', 'normal')
 set(gca, 'XScale', 'log')
title({'Sputtering Yield D on Al','As a Function of Energy and Angle'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label
set(gca,'fontsize',16)
br_edge = interpn(r_grid,z_grid,br_interp',0*z_grid + radius, z_grid);
bz_edge = interpn(r_grid,z_grid,bz_interp',0*z_grid + radius, z_grid);
bfield_angle_z = acosd(br_edge./sqrt(br_edge.^2 + bz_edge.^2));
element_angle = interpn(z_grid,bfield_angle_z,centroid(:,3));
element_angle(find(element_angle >90)) = 180 - element_angle(find(element_angle >90));
element_angle(find(element_angle >=89.9)) = 89.9;
Y0=interpn(energy,angle,spyld',e_value,element_angle','pchip',0);


plotSet = 1:1:length(planes);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
figure(102)
patch(transpose(X),transpose(Y),transpose(Z),element_angle,'FaceAlpha',1,'EdgeColor','none')
colorbar
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

ion_temp_wall = interpn(y,z,ion_temp,0.06256,0);
ion_dens_wall = interpn(y,z,dens,0.06256,0);
ion_dens_wall = dens_value;

k=1.38e-23*11604;
c_bar = sqrt(8*k*ion_temp_wall/pi/4/1.66e-27);
flux = 0.25*ion_dens_wall*c_bar;

plotSet = find(surfs>0);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
% THETA = atan2(Y,X);
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

figure(103)
patch(transpose(THETA),transpose(Z),0*transpose(Z),flux(plotSet).*Y0(plotSet),'FaceAlpha',1,'EdgeColor','k')
% patch(transpose(X),transpose(Y),transpose(Z),flux(plotSet).*Y0(plotSet),'FaceAlpha',1,'EdgeColor','k')

colorbar
title({'D eroded Al Flux [m^{-2}s^{-1}]','Low Flux Case'})
xlabel('Theta [radian]') % x-axis label
ylabel('z [m]') % y-axis label
set(gca,'fontsize',16)

xx = -0.06;
yy = 0.0;
zz = 0.0;
distance = sqrt((centroid(:,1) - xx).^2 + (centroid(:,2) - yy).^2 + (centroid(:,3) - zz).^2);
[v i] = min(distance);

nP = 1000000;
erosion = flux.*Y0.*area';

erosion_inds = find(erosion);
erosion_sub = erosion(erosion_inds);
erosion_sub_cdf = cumsum(erosion_sub);
erosion_rate=erosion_sub_cdf(end);
erosion_sub_cdf = erosion_sub_cdf./erosion_sub_cdf(end);

% plot(erosion_sub_cdf)

rand1 = rand(nP,1);

element = interp1([0, erosion_sub_cdf],0:1:length(erosion_sub_cdf),rand1);

element_ceil = ceil(element);
x_sample = zeros(1,nP);
y_sample = zeros(1,nP);
z_sample = zeros(1,nP);
vx_sample = zeros(1,nP);
vy_sample = zeros(1,nP);
vz_sample = zeros(1,nP);
m = 27;
% meanE = 4.0;
% v = sqrt(2*meanE*1.602e-19/m/1.66e-27);
nPoints = 200;
maxE = 20;
Eb = 3.39;%8.79;
a = 5;
E = linspace(0,maxE,nPoints);
dE = E(2);
thompson2 = a*(a-1)*E.*Eb^(a-1)./(E+Eb).^(a+1);
figure(234)
plot(E,thompson2)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Sample Particle Positions')

m=27;
ecdf = cumsum(thompson2);
ecdf = ecdf./ecdf(end);
rand1 = rand(1,nP);
randTheta = 2*pi*rand(1,nP);
randPhi = 0.5*pi*rand(1,nP);
Esamp = interp1(ecdf,E,rand1);
v = sqrt(2*Esamp*1.602e-19/m/1.66e-27)';
    vx = v'.*cos(randTheta).*sin(randPhi);
    vy = v'.*sin(randTheta).*sin(randPhi);
    vz = v'.*cos(randPhi);
buffer = 1e-5;
plotSet = 1:length(planes);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
particle_counts = histcounts(erosion_inds(element_ceil),0.5:1:(length(plane_norm)+0.5));
nP0 = 0;

for i=1:length(particle_counts)
    if particle_counts(i) > 0
% i
% ind = erosion_inds(element_ceil);
 x_tri = X(i,:);
 y_tri = Y(i,:);
 z_tri = Z(i,:);
    parVec = [x_tri(2) - x_tri(1), y_tri(2) - y_tri(1) , z_tri(2) - z_tri(1)];
    parVec = parVec./norm(parVec);
samples = sample_triangle(x_tri,y_tri,z_tri,particle_counts(i));

normal = inDir(i)*(-abcd(i,1:3)/plane_norm(i));

v_inds = nP0+1:nP0+particle_counts(i);

x_sample(v_inds) = samples(:,1) + buffer*normal(1);
y_sample(v_inds) = samples(:,2) + buffer*normal(2);
z_sample(v_inds) = samples(:,3) + buffer*normal(3);

parVec2 = cross(parVec,normal);

newV = vx(v_inds)'.*parVec + vy(v_inds)'.*parVec2 + vz(v_inds)'.*normal;
% newV = v(i)*normal;
vx_sample(v_inds) = newV(:,1);
vy_sample(v_inds) = newV(:,2);
vz_sample(v_inds) = newV(:,3);

nP0 = nP0+particle_counts(i)
    end
end

index_array = 1:1:nP;
index_array(randperm(length(index_array)));
x_sample = x_sample(index_array);
y_sample = y_sample(index_array);
z_sample = z_sample(index_array);
vx_sample = vx_sample(index_array);
vy_sample = vy_sample(index_array);
vz_sample = vz_sample(index_array);

figure(103)
hold on
theta_sample = atan2(y_sample,x_sample);
% scatter3(x_sample,y_sample,z_sample,'k')
% scatter3(theta_sample,z_sample,0*z_sample,'r')
% quiver3(x_sample,y_sample,z_sample,vx_sample./10./v,vy_sample./10./v,vz_sample./10./v)
% quiver3(0*vx',0*vx',0*vx',vx'./v,vy'./v,vz'./v)

xlabel('Theta [degrees]')
ylabel('z [m]')
zlabel('Z [m]')
title('Eroded Flux')
axis([0 360 -0.2 0.2])

ncid = netcdf.create(['./particle_source_helicon.nc'],'NC_WRITE')
 
dimP = netcdf.defDim(ncid,'nP',nP);

xVar = netcdf.defVar(ncid,'x','double',[dimP]);
yVar = netcdf.defVar(ncid,'y','double',[dimP]);
zVar = netcdf.defVar(ncid,'z','double',[dimP]);
vxVar = netcdf.defVar(ncid,'vx','double',[dimP]);
vyVar = netcdf.defVar(ncid,'vy','double',[dimP]);
vzVar = netcdf.defVar(ncid,'vz','double',[dimP]);

netcdf.endDef(ncid);
 
netcdf.putVar(ncid, xVar, x_sample);
netcdf.putVar(ncid, yVar, y_sample);
netcdf.putVar(ncid, zVar, z_sample);
netcdf.putVar(ncid, vxVar, vx_sample);
netcdf.putVar(ncid, vyVar, vy_sample);
netcdf.putVar(ncid, vzVar, vz_sample);

netcdf.close(ncid);
%% check erosion pattern
% theta_sample = atan2d(y_sample,x_sample)+180;
% 
% walls = find(surfs ==1);
% 
% wall_bins = zeros(1,length(walls));
% for i=1:length(wall_bins)
%     ind = walls(i);
%   in = inpolygon(theta_sample,z_sample,THETA(i,:),Z(ind,:));
%   wall_bins(i) = sum(length(find(in)));
% end
% figure(201)
% patch(transpose(THETA),transpose(Z(walls,:)),0*transpose(Z(walls,:)),wall_bins,'FaceAlpha',1,'EdgeColor','k')
% % patch(transpose(X),transpose(Y),transpose(Z),flux(plotSet).*Y0(plotSet),'FaceAlpha',1,'EdgeColor','k')
% 
% colorbar
% title({'Sampled Erosion'})
% xlabel('Theta [radian]') % x-axis label
% ylabel('z [m]') % y-axis label
% set(gca,'fontsize',16)
% save('erosion0.mat','wall_bins');
function samples = sample_triangle(x,y,z,nP)
x_transform = x - x(1);
y_transform = y - y(1);
z_transform = z- z(1);

% figure(2)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])

v1 = [x_transform(2) y_transform(2) z_transform(2)];
v2 = [x_transform(3) y_transform(3) z_transform(3)];
v12 = v2 - v1;
normalVec = cross(v1,v2);

a1 = rand(nP,1);
a2 = rand(nP,1);

samples = a1.*v1 + a2.*v2;
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
samples2x = samples(:,1) - v2(1);
samples2y = samples(:,2) - v2(2);
samples2z = samples(:,3) - v2(3);
samples12x = samples(:,1) - v1(1);
samples12y = samples(:,2) - v1(2);
samples12z = samples(:,3) - v1(3);
v1Cross = [(v1(2).*samples(:,3) - v1(3).*samples(:,2)) (v1(3).*samples(:,1) - v1(1).*samples(:,3)) (v1(1).*samples(:,2) - v1(2).*samples(:,1))];
v2 = -v2;
v2Cross = [(v2(2).*samples2z - v2(3).*samples2y) (v2(3).*samples2x - v2(1).*samples2z) (v2(1).*samples2y - v2(2).*samples2x)];
v12Cross = [(v12(2).*samples12z - v12(3).*samples12y) (v12(3).*samples12x - v12(1).*samples12z) (v12(1).*samples12y - v12(2).*samples12x)];

v1CD = normalVec(1)*v1Cross(:,1) + normalVec(2)*v1Cross(:,2) + normalVec(3)*v1Cross(:,3);
v2CD = normalVec(1)*v2Cross(:,1) + normalVec(2)*v2Cross(:,2) + normalVec(3)*v2Cross(:,3);
v12CD = normalVec(1)*v12Cross(:,1) + normalVec(2)*v12Cross(:,2) + normalVec(3)*v12Cross(:,3);

inside = abs(sign(v1CD) + sign(v2CD) + sign(v12CD));
insideInd = find(inside ==3);
notInsideInd = find(inside ~=3);
% scatter3(samples(insideInd,1),samples(insideInd,2),samples(insideInd,3))

 v2 = -v2;
dAlongV1 = v1(1).*samples(notInsideInd,1) + v1(2).*samples(notInsideInd,2) + v1(3).*samples(notInsideInd,3);
dAlongV2 = v2(1).*samples(notInsideInd,1) + v2(2).*samples(notInsideInd,2) + v2(3).*samples(notInsideInd,3);

dV1 = norm(v1);
dV2 = norm(v2);
halfdV1 = 0.5*dV1;
halfdV2 = 0.5*dV2;

samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v1(1))+0.5*v1(1) ...
    -(samples(notInsideInd,2) - 0.5*v1(2))+0.5*v1(2) ...
    -(samples(notInsideInd,3) - 0.5*v1(3))+0.5*v1(3)];
% samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v2(1))+0.5*v2(1) ...
%     -(samples(notInsideInd,2) - 0.5*v2(2))+0.5*v2(2) ...
%     -(samples(notInsideInd,3) - 0.5*v2(3))+0.5*v2(3)];
samples(notInsideInd,:) = [(samples(notInsideInd,1) + v2(1)) ...
    (samples(notInsideInd,2) +v2(2)) ...
    (samples(notInsideInd,3) + v2(3))];
% figure(4)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))

samples(:,1) = samples(:,1)+ x(1);
samples(:,2) = samples(:,2)+ y(1);
samples(:,3) = samples(:,3)+ z(1);

% figure(5)
% plot3([x x(1)],[y y(1)],[z z(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
end
