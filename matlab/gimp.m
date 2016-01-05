clear variables

% Load physical constants
% SI Units throughout

constants

% Load run input file

init2

% Load ADAS data (cross sections, ionization rates, etc)

[Tb_inz, dens_inz, RateCoeff_inz, State_inz] = ADF11(file_inz);
[Tb_rcmb, dens_rcmb, RateCoeff_rcmb, State_rcmb] = ADF11(file_rcmb);

% Create volume grid

yMinV = yMin;
yMaxV = yMax;

zMinV = zMin;
zMaxV = zMax;

xV_1D = linspace(xMinV,xMaxV,nXv);
yV_1D = linspace(yMinV,yMaxV,nYv);
zV_1D = linspace(zMinV,zMaxV,nZv);

dXv = xV_1D(2)-xV_1D(1);
dYv = yV_1D(2)-yV_1D(1);
dZv = zV_1D(2)-zV_1D(1);

global xyz 
xyz.x = xV_1D; % Create volume coordinate strucutre
xyz.y = yV_1D;
xyz.z = zV_1D;

% Create surface grid

surf_y1D = linspace(yMin,yMax,nY);
surf_z1D = linspace(zMin,zMax,nZ);
surf_x2D = zeros(nY,nZ);
surf_hist = zeros(nY,nZ);

for j=1:nY 
    surf_x2D(j,:) =  (surf_z1D - surface_zIntercept) / surface_dz_dx;
end

% Plot initial (zeroed) surface histogram

if plotInitialSurface 
    figure(1)
   h1 =  surf(surf_z1D,surf_y1D,surf_x2D,surf_hist);

            %xlabel('z axis')
            %ylabel('y axis')
            %zlabel('x axis')
            title('Surface')
            axis equal
            drawnow
            xlim([xMinV xMaxV])
            ylim([yMinV yMaxV])
            zlim([zMinV zMaxV])
            az = 90;
            el = 90;
            view(az,el);
            %set(gca,'Zdir','reverse')
end



% Setup B field

Bx = zeros(nXv,nYv,nZv);
By = zeros(nXv,nYv,nZv);
Bz = zeros(nXv,nYv,nZv);

Bx(:) = Bx_in;
By(:) = By_in;
Bz(:) = Bz_in;

Bfield.x = Bx; % Create background B strucutre
Bfield.y = By;
Bfield.z = Bz;

BMag = sqrt( Bx.^2 + By.^2 + Bz.^2 );

% Setup profiles

[n_ nS] = size(background_amu);

density_m3 = zeros(nXv,nYv,nZv,nS);
temp_eV = zeros(nXv,nYv,nZv,nS);

Ex = zeros(nXv,nYv,nZv); 
Ey = zeros(nXv,nYv,nZv);
Ez = zeros(nXv,nYv,nZv);

Efield_3D.x = Ex; % Create E field structure
Efield_3D.y = Ey;
Efield_3D.z = Ez;


perDiffusionCoeff = zeros(nXv,nYv,nZv);

V_1D = sheathPotential * exp( xV_1D / sheathWidth );

for i=1:nXv
    for j=1:nYv
        for k=1:nZv
            
            density_m3(i,j,k,:) = maxDensity * exp( (xMinV-xV_1D(i)) / densitySOLDecayLength );
            temp_eV(i,j,k,:) = maxTemp_eV * exp( (xMinV-xV_1D(i)) / tempSOLDecayLength );
            
            perDiffusionCoeff(i,j,k) = perDiffusionCoeff_in;
            
            if i>1 && i<nXv
                Efield_3D.x(i,j,k) = -(V_1D(i+1)-V_1D(i-1)) / (2*dXv);
            elseif i==1
                Efield_3D.x(i,j,k) = -(-1*V_1D(i)+V_1D(i+1)) / dXv;
            elseif i==nXv
                Efield_3D.x(i,j,k) = -(-V_1D(i-1)+V_1D(i)) / dXv;
            end
            
            Efield_3D.y(i,j,k) = 0;
            Efield_3D.z(i,j,k) = 0;
            
        end
    end
end

% Plot slices through the profiles
if plot1DProfileSlices
    figure(2)
    subplot(2,2,1)
    semilogy(xV_1D,density_m3(:,1,1,1))
    title('Electron density')
    subplot(2,2,2)
    semilogy(xV_1D,temp_eV(:,1,1,1))
    title('Electron temp [eV]')
    subplot(2,2,3)
    plot(xV_1D,V_1D(:))
    title('Sheath potential [V]')
    subplot(2,2,4)
    plot(xV_1D,Ex(:,1,1))
    title('Ex [V/m]')
end

% Populate the impurity particle list

for p=1:nP
   particles(p) = particle;
   
   particles(p).Z = impurity_Z;
   particles(p).amu = impurity_amu;
   
   particles(p).x = x_start;
   particles(p).y = y_start;
   particles(p).z = z_start;
   
   particles(p).vx = sign(energy_eV_x_start) * sqrt(2*abs(energy_eV_x_start*Q)/(particles(p).amu*MI));
   particles(p).vy = sign(energy_eV_y_start) * sqrt(2*abs(energy_eV_y_start*Q)/(particles(p).amu*MI));
   particles(p).vz = sign(energy_eV_z_start) * sqrt(2*abs(energy_eV_z_start*Q)/(particles(p).amu*MI));
  
   particles(p).hitWall = 0;
   
   [s1,s2,s3,s4,s5,s6] = RandStream.create('mrg32k3a','NumStreams',6,'Seed','shuffle'); %Include ,'Seed','shuffle' to get different values each time

   particles(p).streams.ionization = s1;
   particles(p).streams.recombination = s2;
   particles(p).streams.perDiffusion = s3;
   particles(p).streams.parVelocityDiffusion = s4;
   particles(p).streams.per1VelocityDiffusion = s5;
   particles(p).streams.per2VelocityDiffusion = s6;

end

% Calculate time step (dt)

max_B = max( BMag(:) );
min_m = impurity_amu * MI;
max_wc = impurity_Z * Q * max_B / min_m;
dt = 2 * pi / max_wc / nPtsPerGyroOrbit;

% Setup arrays to store history

pre_history    

% Main loop

t = cputime;
parfor p=1:nP
    
    for n_steps = 1:nT
            
        time = n_steps * dt;
         
        [T_local, n_local] = particles(p).Tn_interp(xyz,temp_eV,density_m3,nS);

        if mod(n_steps, ionization_nDtPerApply) == 0
            n_ionization_apply = n_steps/ionization_nDtPerApply;
            particles(p).ionization(ionization_nDtPerApply*dt,T_local,n_local,RateCoeff_inz,Tb_inz, dens_inz,State_inz,rand(particles(p).streams.ionization));
            particles(p).recombination(ionization_nDtPerApply*dt,T_local,n_local,RateCoeff_rcmb,Tb_rcmb,dens_rcmb,State_rcmb,rand(particles(p).streams.recombination));
        end
      
        [nu_s, nu_d, nu_par,nu_E] = particles(p).slow(T_local,n_local,nS,background_amu,background_Z);
  
        [E_local, B_local] = particles(p).field_interp(xyz,Bfield,Efield_3D);

        [e1, e2, e3] = particles(p).direction(B_local,E_local);
        
        particles(p).cfDiffusion(B_local,xyz,perDiffusionCoeff,dt,rand(particles(p).streams.perDiffusion));
        
        diagnostics = particles(p).dv_coll(e1,e2,e3,nu_s,nu_d,nu_par,nu_E,dt,rand(particles(p).streams.parVelocityDiffusion), rand(particles(p).streams.per1VelocityDiffusion), rand(particles(p).streams.per2VelocityDiffusion));
        
        % Boris integrator
        
        particles(p).boris(xyz,Efield_3D,B_local,dt);

        %history
        xHistory(n_steps,p) = particles(p).x;
        yHistory(n_steps,p) = particles(p).y;
        zHistory(n_steps,p) = particles(p).z;
       

        % find particles that returned to the wall

        if particles(p).x > ( particles(p).z - surface_zIntercept ) / surface_dz_dx;
           
            particles(p).amu = 1e20;
            particles(p).vx = 0;
            particles(p).vy = 0;
            particles(p).vz = 0;
            particles(p).hitWall = 1;
            
        end

    end
end

history_plot

e = cputime-t

surface_scatter

