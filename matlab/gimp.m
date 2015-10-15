clear variables

% SI Units throughout

constants

% Surface grid

nY = 9;%Number of surface cells
nZ = 10;

yMin = -0.05;
yMax = +0.05;

zMin = -0.05;
zMax = +0.05;

y = linspace(yMin,yMax,nY);
z = linspace(zMin,zMax,nZ);

% Volume grid

nXv = 100;
nYv = 80;
nZv = 50;

xMinV = -0.1;
xMaxV = 0;

yMinV = yMin;
yMaxV = yMax;

zMinV = zMin;
zMaxV = zMax;

xV = linspace(xMinV,xMaxV,nXv);
yV = linspace(yMinV,yMaxV,nYv);
zV = linspace(zMinV,zMaxV,nZv);

dXv = xV(2)-xV(1);
dYv = yV(2)-yV(1);
dZv = zV(2)-zV(1);

% Setup B field

Bx = zeros(nXv,nYv,nZv);
By = zeros(nXv,nYv,nZv);
Bz = zeros(nXv,nYv,nZv);

Bx(:) = 0.0;
By(:) = 1.0;
Bz(:) = 0.14;

Bfield.x = Bx;
Bfield.y = By;
Bfield.z = Bz;

BMag = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Setup profiles

Z = [-1 1]
amu = [ME/MI 1];

[n_ nS] = size(amu);

density = zeros(nXv,nYv,nZv,nS);
temp_eV = zeros(nXv,nYv,nZv,nS);

Ex = zeros(nXv,nYv,nZv); 
Ey = zeros(nXv,nYv,nZv);
Ez = zeros(nXv,nYv,nZv);

Efield.x = Ex;
Efield.y = Ey;
Efield.z = Ez;

xyz.x = xV;
xyz.y = yV;
xyz.z = zV;

maxDensity = 1e19;
densitySOLDecayLength = 0.01;

maxTemp_eV = 20;
tempSOLDecayLength = 0.02;

sheathPotential = -60;
sheathWidth = 0.005;

V = sheathPotential * exp( xV / sheathWidth );

for i=1:nXv
    for j=1:nYv
        for k=1:nZv
            
            density(i,j,k,:) = maxDensity * exp( (xMinV-xV(i)) / densitySOLDecayLength );
            temp_eV(i,j,k,:) = maxTemp_eV * exp( (xMinV-xV(i)) / tempSOLDecayLength );
            
            if i>1 && i<nXv
                Ex(i,j,k) = (V(i+1)-V(i-1)) / (2*dXv);
            elseif i==1
                Ex(i,j,k) = (-1*V(i)+V(i+1)) / dXv;
            elseif i==nXv
                Ex(i,j,k) = (-V(i-1)+V(i)) / dXv;
            end
            Ey(i,j,k) = 0;
            Ez(i,j,k) = 0;
            
        end
    end
end

subplot(2,2,1)
semilogy(xV,density(:,1,1,1))
title('Electron density')
subplot(2,2,2)
semilogy(xV,temp_eV(:,1,1,1))
title('Electron temp [eV]')
subplot(2,2,3)
plot(xV,V(:))
title('Sheath potential [V]')
subplot(2,2,4)
plot(xV,Ex(:,1,1))
title('Ex [V/m]')


% Populate the particle list

nP = 10000;
yTileMin = -0.005;
yTileMax = +0.005;
zTileMin = -0.005;
zTileMax = +0.005;

bindingEnergy_eV = 8;
maxEnergy_eV = 100;
[energy_x, energy_y, energy_z] = thompson_energyXYZ( bindingEnergy_eV, maxEnergy_eV, nP);

for p=1:nP
   particles(p) = particle;
   
   particles(p).Z = 1;
   particles(p).amu = 184;
   
   particles(p).x = 0;
   particles(p).y = (rand * (yTileMax-yTileMin) ) * yTileMin;
   particles(p).z = (rand * (zTileMax-zTileMin) ) * zTileMin;
   
   particles(p).vx = sign(energy_x(p)) * sqrt(2*abs(energy_x(p)*Q)/(particles(p).amu*MI));
   particles(p).vy = sign(energy_y(p)) * sqrt(2*abs(energy_y(p)*Q)/(particles(p).amu*MI));
   particles(p).vz = sign(energy_z(p)) * sqrt(2*abs(energy_z(p)*Q)/(particles(p).amu*MI));
   
end

% N = Npol*Ntor;%Total number of cells
% A = Lt*Lp;%Area of each cell
% 
% Pp = 20;%Number of particles launced per surface cell
% 
% Tij = zeros(N,N);

%tracker_param

nT = 1e2;

max_Z = 3;
max_B = max( BMag(:) );
min_m = 184 * MI;
max_wc = max_Z*Q * max_B / min_m;

nPtsPerGyroOrbit = 4;
dt = 2*pi/max_wc / nPtsPerGyroOrbit;

for p = 1:nP
    
    [t,y] = particles(p).move(nT*dt,dt,Efield,Bfield,xyz);
    
    plot3(y(:,1),y(:,2),y(:,3))
    
end

colormap('winter')
imagesc(Tij)
colorbar

map2 = transpose(reshape(Tij(37,:), [Npol,Ntor]));

imagesc(map2)
colorbar