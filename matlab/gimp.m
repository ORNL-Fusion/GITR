clear variables

% SI Units throughout

constants
init2


[Tb_inz, dens_inz, RateCoeff_inz, State_inz] = ADF11(file_inz);

[Tb_rcmb, dens_rcmb, RateCoeff_rcmb, State_rcmb] = ADF11(file_rcmb);
clearvars file_inz file_rcmb;

% Surface grid
if exist('nY','var') == 0
    nY = 9;%Number of surface cells
    nZ = 10;
    
    yMin = -0.05;
    yMax = +0.05;
    
    zMin = -0.05;
    zMax = +0.05;
    surf_slope = 0;
    surf_incpt = 0;
end
y = linspace(yMin,yMax,nY);
z = linspace(zMin,zMax,nZ);
surfx = zeros(nY,nZ);
surf_hist = surfx;
for k=1:nZ 
surfx(:,k) = surf_slope*y + surf_incpt;
end
figure(3)
surf(z,y,surfx,surf_hist)
            xlabel('z axis')
            ylabel('y axis')
            zlabel('x axis')
            title('Surface')
            axis equal

% Volume grid
if exist('nXv','var') == 0
    nXv = 100;
    nYv = 80;
    nZv = 50;
    
    xMinV = -0.005;
    xMaxV = 0;
end

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

Bx(:) = -0.2;
By(:) = 1.0;
Bz(:) = 0.3;

if exist('Bfieldx_dat','var')==1
    Bx(:) = Bfieldx_dat;
    By(:) = Bfieldy_dat;
    Bz(:) = Bfieldz_dat;
end

Bfield.x = Bx;
Bfield.y = By;
Bfield.z = Bz;

BMag = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Setup profiles

Z = [-1 1];
amu = [ME/MI 2];

[n_ nS] = size(amu);

density = zeros(nXv,nYv,nZv,nS);
temp_eV = zeros(nXv,nYv,nZv,nS);

Ex = zeros(nXv,nYv,nZv); 
Ey = zeros(nXv,nYv,nZv);
Ez = zeros(nXv,nYv,nZv);

Dperp = zeros(nXv,nYv,nZv);

Efield.x = Ex;
Efield.y = Ey;
Efield.z = Ez;

xyz.x = xV;
xyz.y = yV;
xyz.z = zV;

if exist('maxDensity','var')==0
    maxDensity = 1e19;
    densitySOLDecayLength = .01;
    
    maxTemp_eV = 20;
    tempSOLDecayLength = 1e-2;
    
    sheathPotential = -60;
    sheathWidth = 0.0001;
    
    Dperp_dat = 0.4;
end

V = sheathPotential * exp( xV / sheathWidth );

for i=1:nXv
    for j=1:nYv
        for k=1:nZv
            
            density(i,j,k,:) = maxDensity * exp( (xMinV-xV(i)) / densitySOLDecayLength );
            temp_eV(i,j,k,:) = maxTemp_eV * exp( (xMinV-xV(i)) / tempSOLDecayLength );
            
            Dperp(i,j,k) = Dperp_dat;
            
            if i>1 && i<nXv
                Ex(i,j,k) = -(V(i+1)-V(i-1)) / (2*dXv);
            elseif i==1
                Ex(i,j,k) = -(-1*V(i)+V(i+1)) / dXv;
            elseif i==nXv
                Ex(i,j,k) = -(-V(i-1)+V(i)) / dXv;
            end
            Ey(i,j,k) = 0;
            Ez(i,j,k) = 0;
            
        end
    end
end
figure(1)
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
if exist('nP','var')==0
    nP = 1;
    Ex_dat = 1;
    Ey_dat = -1;
    Ez_dat = -1;
    x_dat = -0.004;
    y_dat = 0;
    z_dat = 0;
    charge = 1;
    mass = 12;
end

yTileMin = -0.005;
yTileMax = +0.005;
zTileMin = -0.005;
zTileMax = +0.005;

%bindingEnergy_eV = 8;
%maxEnergy_eV = 100;
%[energy_x, energy_y, energy_z] = thompson_energyXYZ( bindingEnergy_eV, maxEnergy_eV, nP);
energy_x(1:nP) = Ex_dat;
energy_y(1:nP) = Ey_dat;
energy_z(1:nP) = Ez_dat;
x_start = x_dat;
y_start = y_dat;
z_start = z_dat;
charge = charge_dat;
mass = mass_dat;
for p=1:nP
   particles(p) = particle;
   
   particles(p).Z = charge;
   particles(p).amu = mass;
   
   particles(p).x = x_start;
   particles(p).y = y_start;%(rand(s1) * (yTileMax-yTileMin) ) * yTileMin;
   particles(p).z = z_start;%(rand(s2) * (zTileMax-zTileMin) ) * zTileMin;
   
   particles(p).vx = sign(energy_x(p)) * sqrt(2*abs(energy_x(p)*Q)/(particles(p).amu*MI));
   particles(p).vy = sign(energy_y(p)) * sqrt(2*abs(energy_y(p)*Q)/(particles(p).amu*MI));
   particles(p).vz = sign(energy_z(p)) * sqrt(2*abs(energy_z(p)*Q)/(particles(p).amu*MI));
  
   particles(p).hitWall = 0;
end



if exist('nT','var')==0
    nT = 1e3;
    max_nT = nT;
    
    max_Z = 3;
    max_B = max( BMag(:) );
    min_m = 184 * MI;
    max_wc = max_Z*Q * max_B / min_m;
    
    nPtsPerGyroOrbit = 1e3;
    dt = 2*pi/max_wc / nPtsPerGyroOrbit;
    
    ionization_factor = 1e2;
end

pre_history

set_rand
for n_steps = 1:nT
    time = n_steps*dt;

    for p=1:nP
        
        
        [T_local, n_local] = particles(p).Tn_interp(xyz,temp_eV,density,nS);
        if mod(n_steps, ionization_factor) == 0
            particles(p).ionization(ionization_factor*dt,T_local,n_local,RateCoeff_inz,Tb_inz, dens_inz,State_inz,random(p,n_steps,1));
            particles(p).recombination(ionization_factor*dt,T_local,n_local,RateCoeff_rcmb,Tb_rcmb,dens_rcmb,State_rcmb,random(p,n_steps,2));
        end
      
        [nu_s, nu_d, nu_par,nu_E] = particles(p).slow(T_local,n_local,nS,amu,Z);

        
        [E_local, B_local] = particles(p).field_interp(xyz,Bfield,Efield);

        [e1, e2, e3] = particles(p).direction(B_local,E_local);
        
        particles(p).cfDiffusion(B_local,xyz,Dperp,dt,random(p,n_steps,3));
        
        
 
        diagnostics = particles(p).dv_coll(e1,e2,e3,nu_s,nu_d,nu_par,nu_E,dt,random(p,n_steps,4:6));
        

        
        

        %%%%%%Boris integrator
        
        particles(p).boris(B_local,E_local,dt);

        history
        
        
        
    end
history_plot

    
    % find particles that returned to the wall

particlesWall = [];
for p = 1:nP
    if particles(p).x > surf_slope*particles(p).y + surf_incpt;%-dXv/10
        particlesWall = [particlesWall particles(p)];
        particles(p).amu = 1e20;
        particles(p).vx = 0;
        particles(p).vy = 0;
        particles(p).vz = 0;
        particles(p).hitWall = 1;
    end
end
end


surface_scatter

