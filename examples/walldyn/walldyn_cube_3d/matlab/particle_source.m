close all
clear all

% Read in 3D surface mesh geometry file
if (exist('x1') == 0)
    fid = fopen(strcat(pwd,'/gitrGeometry.cfg'));

    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:18
        tline = fgetl(fid);
        evalc(tline);
    end
    Zsurface = Z;
end
abcd = [a' b' c' d'];
surface = find(Zsurface);
nSurfaces = length(a);

subset = 1:length(x1);

figure
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface)
title('Geometry')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

mass = 174;
index = 1;
nP = 1e4;

        nPoints = 200;
        maxE = 20;
        Eb = 8.79;
        a = 5;
        E = linspace(0,maxE,nPoints);
        dE = E(2);
        thompson = a*(a-1)*E.*Eb^(a-1)./(E+Eb).^(a+1);
        figure
        plot(E,thompson)
        xlabel('E [eV]')
        ylabel('pdf')



        ecdf = cumsum(thompson);
        ecdf = ecdf./ecdf(end);
        rand1 = rand(1,nP);
        randTheta = 2*pi*rand(1,nP);
        randPhi = 0.5*pi*rand(1,nP);
        Esamp = interp1(ecdf,E,rand1);
        v = sqrt(2*Esamp*1.602e-19/mass/1.66e-27)';
        vx = v'.*cos(randTheta).*sin(randPhi);
        vy = v'.*sin(randTheta).*sin(randPhi);
        vz = v'.*cos(randPhi);



        buffer = 1e-6;
        plotSet = 1:length(area);
        planes = [x1' y1' z1' x2' y2' z2' x3' y3' z3'];

        X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
        Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
        Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

        % particle_counts = histcounts(erosion_inds(element_ceil),0.5:1:(length(plane_norm)+0.5));
        nP0 = 0;
        nTriangles = length(planes);

        r1 = rand(nP,1);
        this_triangle = index;
        for j=1:nP

            i = this_triangle;
            x_tri = X(i,:);
            y_tri = Y(i,:);
            z_tri = Z(i,:);
            parVec = [x_tri(2) - x_tri(1), y_tri(2) - y_tri(1) , z_tri(2) - z_tri(1)];
            parVec = parVec./norm(parVec);
            samples = sample_triangle(x_tri,y_tri,z_tri,1);

            normal = -inDir(i)*(-abcd(i,1:3)./plane_norm(i));

            v_inds =j; % nP0+1:nP0+nP;

            x_sample(v_inds) = samples(:,1) + buffer*normal(1);
            y_sample(v_inds) = samples(:,2) + buffer*normal(2);
            z_sample(v_inds) = samples(:,3) + buffer*normal(3);

            parVec2 = cross(parVec,normal);

            newV = vx(v_inds)'.*parVec + vy(v_inds)'.*parVec2 + vz(v_inds)'.*normal;
            % newV = v(i)*normal;
            vx_sample(v_inds) = newV(:,1);
            vy_sample(v_inds) = newV(:,2);
            vz_sample(v_inds) = newV(:,3);


        end
        figure(1)
        hold on
        scatter3(x_sample,y_sample,z_sample)
        quiver3(x_sample,y_sample,z_sample,vx_sample./10./v',vy_sample./10./v',vz_sample./10./v')
        ncid = netcdf.create('particle_source.nc','NC_WRITE')

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

axis equal

function samples = sample_triangle(x,y,z,nP)
x_transform = x - x(1);
y_transform = y - y(1);
z_transform = z - z(1);

v1 = [x_transform(2) y_transform(2) z_transform(2)];
v2 = [x_transform(3) y_transform(3) z_transform(3)];
v12 = v2 - v1;
normalVec = cross(v1,v2);

a1 = rand(nP,1);
a2 = rand(nP,1);

samples = a1.*v1 + a2.*v2;

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

samples(notInsideInd,:) = [(samples(notInsideInd,1) + v2(1)) ...
    (samples(notInsideInd,2) + v2(2)) ...
    (samples(notInsideInd,3) + v2(3))];


samples(:,1) = samples(:,1)+ x(1);
samples(:,2) = samples(:,2)+ y(1);
samples(:,3) = samples(:,3)+ z(1);

end
