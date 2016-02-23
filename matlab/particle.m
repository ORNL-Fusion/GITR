classdef particle < handle
    
    properties
        x
        y
        z
        vx
        vy
        vz
        xPrevious
        yPrevious
        zPrevious
        vxPrevious
        vyPrevious
        vzPrevious
        energy
        Z
        amu
        hitWall
        leftVolume
        streams
        perpDistanceToSurface
        impactAngle
    end
    
    methods
        
        function [T Y] = move (this, end_t, dt, E, B, xyz,interpolators,...
                xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV,...
                surface_zIntercept,surface_dz_dx, ...
                decayLength, potential,background_Z,background_amu,maxTemp_eV)
            
            
            status = 0;
            
            IC = [this.x this.y this.z this.vx this.vy this.vz]';
            tSpan = [0,end_t];
            
            maxStep = dt;
            if this.perpDistanceToSurface < 1e-3 && this.hitWall == 0 && this.leftVolume ==0
                maxStep = 5e-12;
            end
            eventHandle = @(t,y)myEvent(t,y,surface_dz_dx,surface_zIntercept);
            options = odeset('Events',eventHandle,'InitialStep',maxStep,'MaxStep',maxStep); 
            %options = odeset('RelTol', 1e-2);
            Einterpolator = interpolators{1};
            Binterpolator = interpolators{2};
            [T ,Y, TE, YE, IE] = ode45(@(t,y) myode(t,y,this,E,B,xyz,Einterpolator,Binterpolator, ...
                 decayLength, potential,surface_dz_dx, ...
               background_Z,background_amu,maxTemp_eV),tSpan,IC,options);
            if TE
                stop
            end
            
%             n_steps = length(T);
%             index = 0;
%             while this.hitWall == 0 && this.leftVolume ==0 && index < n_steps
%                 index = index+1;
%                 this.x = Y(index,1); % this doesn't seem to work, i.e., cannot modify p outside this scope :(
%                 this.y = Y(index,2);
%                 this.z = Y(index,3);
%                 
%                 this.vx = Y(index,4);
%                 this.vy = Y(index,5);
%                 this.vz = Y(index,6);
%                 this.OutOfDomainCheck(xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV);
%                 this.HitWallCheck(surface_zIntercept,surface_dz_dx);
%                 
%             end
%             if index ~= n_steps
%                 stop
%             end
            this.x = Y(end,1); % this doesn't seem to work, i.e., cannot modify p outside this scope :(
            this.y = Y(end,2);
            this.z = Y(end,3);
            
            this.vx = Y(end,4);
            this.vy = Y(end,5);
            this.vz = Y(end,6);
            
        end
        
        
        
        function borisMove(this,xyz,Efield3D,Bfield3D,dt,interpolators,...
                positionStepTolerance,dv_threshold, ...
                decayLength, potential,surface_dz_dx,surface_zIntercept,background_Z, ...
                background_amu,maxTemp_eV, sheath_timestep_factor)
            if this.hitWall == 0 && this.leftVolume ==0
            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
            time = 0;
            tspan = dt;
            tolerance = 1e-16;
            
            surfaceDirection = [surface_dz_dx 0 -1];
            surfaceDirection_unit = surfaceDirection/norm(surfaceDirection);
            
            perpDistanceToSurface1 = ( -surface_dz_dx*this.x + this.z + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
            while abs(tspan - time) > tolerance && perpDistanceToSurface1 >= 0

                
                BfieldInterpolator = interpolators{2};
                B = BfieldInterpolator(this,xyz,Bfield3D);
                
                EfieldInterpolator = interpolators{1};
                E = EfieldInterpolator(this,xyz,Efield3D, decayLength, potential,surface_dz_dx,B, ...
                    background_Z,background_amu,maxTemp_eV);
                
                BMagPart =norm(B);
                
                % Constants used in Boris method Lorentz Integrator
                q_prime = this.Z*Q/(this.amu*MI)*dt/2;
                coeff = 2*q_prime/(1+(q_prime*BMagPart).^2);
                
                % Boris Method Lorentz Integrator
                v = [this.vx this.vy this.vz];
                r = [this.x this.y this.z];
                v_start = v;
                v_minus = v + q_prime*E;
                
                v = v_minus + q_prime*[v_minus(2)*B(3) - v_minus(3)*B(2), v_minus(3)*B(1) - v_minus(1)*B(3),v_minus(1)*B(2) - v_minus(2)*B(1)];
                
                v = v_minus + coeff*[v(2)*B(3) - v(3)*B(2), v(3)*B(1) - v(1)*B(3),v(1)*B(2) - v(2)*B(1)];
                
                v = v + q_prime*E;
                
                step = v*dt;
                r = r + step;
                
                if abs(norm(v) - norm(v_start))/norm(v_start) > dv_threshold
                    dt = dt/10;
                else
                    perpDistanceToSurface1 = ( -surface_dz_dx*r(1) + r(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
                    if perpDistanceToSurface1 < 0
                        this.hitWall = 1;
                        t = (-surface_zIntercept + surface_dz_dx*this.x - this.z)/(-surface_dz_dx*(this.x -r(1)) + (this.z -r(3)));
                        this.x = this.x + (this.x - r(1))*t;
                        this.y = this.y + (this.y - r(2))*t;
                        this.z = this.z + (this.z - r(3))*t;
                        this.impactAngle = acosd((surfaceDirection_unit(1)*v(1) + surfaceDirection_unit(2)*v(2)+ surfaceDirection_unit(3)*v(3))/norm(v));
                        this.vx = v(1);%this.vx + (this.vx - v(1))*v_interpFactor;
                        this.vy = v(2);%this.vy + (this.vy - v(2))*v_interpFactor;
                        this.vz = v(3);%this.vz + (this.vz - v(3))*v_interpFactor;
                        this.vxPrevious = this.vx;
                        this.vyPrevious = this.vy;
                        this.vzPrevious = this.vz;
                        
                    else
                        this.x =r(1);
                        this.y =r(2);
                        this.z =r(3);
                        this.vx = v(1);
                        this.vy = v(2);
                        this.vz = v(3);
                    end
                    this.PerpDistanceToSurface(surface_dz_dx,surface_zIntercept);
                    time = time + dt;
                end
            end
            end
            
            if isnan(this.x) == 1
                stop
            end
        end
        
        
        function [T Yr] = rk4(part,xV,yV,zV,Bx,By,Bz,BMag,Ex,Ey,Ez,end_t,dt,steps)
            
            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
            m = MI*part.amu;
            q_m = part.Z*Q/m;
            E = [0 0 0];
            B = [0 0 0];
            
            v = [part.vx part.vy part.vz];
            r = [part.x part.y part.z];
            
            h = end_t/steps;
            T = zeros(steps+1,1);
            Yr = zeros(steps+1,3);
            Yv = zeros(steps+1,3);
            T(1) = 0;
            Yr(1,:) = r;
            Yv(1,:) = v;
            for j=1:steps
                tj = T(j);
                yjr = Yr(j,:);
                yjv = Yv(j,:);
                [E B] = field_interp(yjr,xV,yV,zV,Bx,By,Bz,Ex,Ey,Ez);
                k1r = h*yjv;
                k1v = h*q_m*(E + cross(yjv,B));
                
                [E B] = field_interp(yjr+k1r/2,xV,yV,zV,Bx,By,Bz,Ex,Ey,Ez);
                k2r = h*(yjv+ k1v/2);
                k2v = h*q_m*(E + cross((yjv+ k1v/2),B));
                
                
                [E B] = field_interp(yjr+k2r/2,xV,yV,zV,Bx,By,Bz,Ex,Ey,Ez);
                k3r = h*(yjv+ k2v/2);
                k3v = h*q_m*(E + cross((yjv+ k2v/2),B));
                
                [E B] = field_interp(yjr+k3r,xV,yV,zV,Bx,By,Bz,Ex,Ey,Ez);
                k4r = h*(yjv+ k3v);
                k4v = h*q_m*(E + cross((yjv+ k3v),B));
                
                Yr(j+1,:) = yjr + (k1r + 2*k2r + 2*k3r + k4r)/6;
                Yv(j+1,:) = yjv + (k1v + 2*k2v + 2*k3v + k4v)/6;
                T(j+1) = h*j;
                
            end
            
        end
        
        function  ionization(this,dt,xyz,density_m3,temp_eV,...
                IonizationRateCoeff,IonizationTemp,IonizationDensity,...
                ChargeState,interpolators, ionizationProbabilityTolerance)
            
            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
            
            if this.Z < max(ChargeState(:));
                
                r1 = rand(this.streams.ionization);
                
                minT = min(IonizationTemp(:));
                
                nS = size(density_m3,4);
                
                T_local = zeros(nS);
                n_local = zeros(nS);
                
                temperatureInterpolator = interpolators{4};
                densityInterpolator = interpolators{5};
                
                for s=1:nS
                    T_local(s) = temperatureInterpolator(this,xyz,temp_eV,s);
                    n_local(s) = densityInterpolator(this,xyz,density_m3,s);
                end
                
                T=T_local(1);
                
                if log10(T) > minT
                    n=n_local(1);
                    
                   % isreal(n) == 0 | isreal(T) == 0
                    Coeff = interpn(IonizationDensity,IonizationTemp,IonizationRateCoeff(:,:,this.Z+1),log10(n/1e6),log10(T),'linear',0);
                    
                    if ( isnan(Coeff) )
                        error('Ionization interpolation out of range')
                    end
                    
                    tion = 1/(10^Coeff*n/1E6);
                    P1 = 1-exp(-dt/tion);
                    
                    if P1 > ionizationProbabilityTolerance
                        error('Ionization probability is too large')
                        P1
                    end
                    
                    if r1 <= P1
                        this.Z = this.Z+1;
                        this.perpDistanceToSurface;
                    end
                end
            end
        end
        
        function recombination(this,dt,xyz,density_m3,temp_eV,...
                RecombinationRateCoeff,RecombinationTemp,RecombinationDensity,...
                ChargeState,interpolators,ionizationProbabilityTolerance)
            
            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
            nS = size(density_m3,4);
            
            T_local = zeros(nS);
            n_local = zeros(nS);
            
                temperatureInterpolator = interpolators{4};
                densityInterpolator = interpolators{5};
                
                for s=1:nS
                    T_local(s) = temperatureInterpolator(this,xyz,temp_eV,s);
                    n_local(s) = densityInterpolator(this,xyz,density_m3,s);
                end
            
            minT = min(RecombinationTemp);
            minN = min(RecombinationDensity);
            T=T_local(1);
            n=n_local(1);
            T = log10(T);
            n = log10(n/1e6);
            if (T > minT) && (n > minN) && this.Z >0
                
                r2 = rand(this.streams.recombination);
                
                Coeff = interpn(RecombinationDensity,RecombinationTemp,RecombinationRateCoeff(:,:,this.Z),n,T);
                
                tion = 1/(10^Coeff*10^n);
                P1 = 1-exp(-dt/tion);
                
                if P1 > ionizationProbabilityTolerance
                    error('Recombination probability is too large')
                    P1
                end
                
                if r2 <= P1
                    this.Z = this.Z-1;
                end
            end
        end
        
        function [nu_s, nu_d, nu_par, nu_E] = slow(this,xyz,density_m3,temp_eV,amu,Z,interpolators,Bfield, ...
                connectionLength,surface_dz_dx,surface_zIntercept,background_amu,flowVelocity_ms)
            
            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
            %persistent v_norm_persistent;
            
            nS = size(density_m3,4);
            
            T_local = zeros(nS);
            n_local = zeros(nS);
            temperatureInterpolator = interpolators{4};
            densityInterpolator = interpolators{5};
            
            for s=1:nS
                T_local(s) = temperatureInterpolator(this,xyz,temp_eV,s);
                n_local(s) = densityInterpolator(this,xyz,density_m3,s);
            end
            
            v = [this.vx this.vy this.vz];
           
            BfieldInterpolator = interpolators{2};
            B_local = BfieldInterpolator(this,xyz,Bfield);
            
            
            for s=1:nS
            flowVelocityInterpolator = interpolators{3};
            flowVelocity = flowVelocityInterpolator(this,xyz,flowVelocity_ms,s,B_local, ...
                connectionLength,temp_eV,background_amu,surface_dz_dx,surface_zIntercept);
            end 
            
            v_relative = v - flowVelocity;
            v_norm = norm(v_relative);
            %             if v_norm == 0
            %                 v_norm = v_norm_persistent;
            %             else
            %                 v_norm_persistent = v_norm;
            %             end
            nu_s = 0;
            nu_d = 0;
            nu_par = 0;
            nu_E = 0;
            
            z = this.Z;
            m = this.amu*MI;
            for j=2:nS
                
                zbackground = Z(j);
                mbackground = amu(j)*MI;
                
                T = T_local(j);
                n = n_local(j);
                
                lam_d = sqrt(EPS0*T/(n*zbackground^2*Q));%only one q in order to convert to J
                lam = 4*pi*n*lam_d^3;
                gam = Q^4*z^2*zbackground^2*log(lam)/(m*m*4*pi*EPS0*EPS0);
                
                a = mbackground/(2*T*Q); %q is just to convert units - no z needed
                
                x = v_norm^2*a;
                psi_prime = 2*sqrt(x/pi)*exp(-x);
                psi_psiprime = erf(sqrt(x));
                psi = psi_psiprime - psi_prime;
                nu_0 = gam*n/v_norm^3;
                nu_s =nu_s -(1+m/mbackground)*psi*nu_0;
                nu_d = nu_d + 2*(psi_psiprime - psi/(2*x))*nu_0;
                nu_par = nu_par + psi/x*nu_0;
                nu_E = nu_E+2*(m/mbackground*psi - psi_prime)*nu_0;
                
            end
            
            
            
        end
        function [e1, e2, e3] = direction(this,flowVelocity,B_local)
            
            
            v = [this.vx this.vy this.vz];
            B_unit = B_local/norm(B_local);
            
            if this.Z ==0
                flowVelocity = [0 0 0];
            end
            
            v_relative = v- flowVelocity ;
            
            g = v_relative;
            e3 = g/norm(g);
            
            s1 = e3(1)*B_unit(1)+e3(2)*B_unit(2)+e3(3)*B_unit(3);
            s2 = sqrt(1-s1*s1);
            
            e1 = 1/s2*(s1*e3 - B_unit);
            e2 = -1/s2*cross(e3,B_unit);
            
            if s2 == 0
                e2 = [e3(3) e3(1) e3(2)];
                s1 = e3(1)*e2(1)+e3(2)*e2(2)+e3(3)*e2(3);
                s2 = sqrt(1-s1*s1);
                e1 = -1/s2*cross(e3,e2);
            end
            
        end

        function CrossFieldDiffusion(this,xyz,Bfield,Dperp,interpolators,dt,positionStepTolerance)
            
            if this.hitWall == 0 && this.leftVolume == 0 && this.Z > 0
                
                s3 = rand(this.streams.perDiffusion);
                
                BfieldInterpolatorHandle = interpolators{2};
                B_local = BfieldInterpolatorHandle(this,xyz,Bfield);
                B_unit = B_local/norm(B_local);
                
                % Calculation of direction perpendicular to B
                
                phi_rnd = 2*pi*s3;
                eperp(1) = cos(phi_rnd);
                eperp(2) = sin(phi_rnd);
                eperp(3) = (-eperp(1)*B_unit(1) - eperp(2)*B_unit(2))/B_unit(3);
                if B_unit(3) == 0
                    eperp(3) = eperp(2);
                    eperp(2) = (-eperp(1)*B_unit(1) - eperp(3)*B_unit(3))/B_unit(2);
                end
                if B_unit == [1 0 0]
                    eperp(3) = eperp(1);
                    eperp(1) = 0;
                elseif B_unit == [0 1 0]
                    eperp(2) = 0;
                elseif B_unit == [0 0 1]
                    eperp(3) = 0;
                end
                norme = norm([eperp(1) eperp(2) eperp(3)]);
                eperp(1) = eperp(1)/norme;
                eperp(2) = eperp(2)/norme;
                eperp(3) = eperp(3)/norme;
                
                DiffusionInterpolatorHandle = interpolators{6};
                D_local = DiffusionInterpolatorHandle(this,xyz,Dperp);
                 step = sqrt(6*D_local*dt);
                if step > positionStepTolerance
                    step
                    error('Position step is too large in cross-field diffusion')
                end
                
                this.x = this.x + step*eperp(1);
                this.y = this.y + step*eperp(2);
                this.z = this.z + step*eperp(3);
            end
            
        end

        function diagnostics = CoulombCollisions(this,xyz,Bfield,flowVelocity_ms,density_m3,temp_eV,...
                background_amu,background_Z,interpolators,...
                dt,velocityChangeTolerance, connectionLength, ...
                surface_dz_dx,surface_zIntercept)
            
            global ME MI Q EPS0;
            
            rand_parVelocityDiffusion = rand(this.streams.parVelocityDiffusion);
            rand_per1VelocityDiffusion = rand(this.streams.per1VelocityDiffusion);
            rand_per2VelocityDiffusion = rand(this.streams.per2VelocityDiffusion);
            
            BfieldInterpolator = interpolators{2};
            B_local = BfieldInterpolator(this,xyz,Bfield);
            
            nS = size(density_m3,4);
            flowVelocityInterpolator = interpolators{3};
            for s=1:nS
                
                flowVelocity = flowVelocityInterpolator(this,xyz,flowVelocity_ms,s,B_local, ...
                    connectionLength,temp_eV,background_amu,surface_dz_dx,surface_zIntercept);
            end
            
            [e1, e2, e3] = this.direction(flowVelocity,B_local);
            
            [nu_s, nu_d, nu_par, nu_E] = this.slow(xyz,density_m3,temp_eV,...
                background_amu,background_Z,interpolators,Bfield,...
                connectionLength,surface_dz_dx,surface_zIntercept,background_amu,flowVelocity_ms);
            v = [this.vx this.vy this.vz];
            v_norm = norm(v);
            
            
            if this.Z ==0
                flowVelocity = [0 0 0];
            end
            v_relative = v- flowVelocity ;
            vRelative_norm = norm(v_relative);
            
            T = v_norm^2*this.amu*MI*pi/8/Q;
            
            dv_slow = e1*(nu_s*dt);
            norm_slow = norm(v_norm*dv_slow);
            
            plus_minus1 = round(rand_parVelocityDiffusion)*2-1;
            dv_par = e1*plus_minus1*sqrt(nu_par*dt);
            dv_parallel = dv_slow+ dv_par;
            norm_par = norm(v_norm*dv_par);
            norm_parallel = norm(v_norm*dv_parallel);
            
            
            plus_minus2 = round(rand_per1VelocityDiffusion)*2-1;
            dv_perp1 =  e1*plus_minus2*sqrt(nu_d/2*dt);
            norm_perp1 = plus_minus2*norm(v_norm*dv_perp1);
            plus_minus3 = round(rand_per2VelocityDiffusion)*2-1;
            dv_perp2 =  e2*plus_minus3*sqrt(nu_d/2*dt);
            norm_perp2 = plus_minus3*norm(v_norm*dv_perp2);
            
            
             ez = (1+nu_s*dt+ plus_minus1*sqrt(nu_par*dt/2) );%
             v_collisions = vRelative_norm*(1-nu_E*dt)*(e3*ez + dv_perp1 + dv_perp2);


            this.vx = v_collisions(1)+ flowVelocity(1);
            this.vy = v_collisions(2)+ flowVelocity(2);
            this.vz = v_collisions(3)+ flowVelocity(3);
            
            dv_collisions = v_collisions - v_norm;

            if dv_collisions > velocityChangeTolerance
                error('Velocity change is too large in Coulomb Collisions')
                dv_collisions
            end

            diagnostics = [T dv_collisions norm_slow norm_par norm_parallel norm_perp1 norm_perp2];
        end
        
        function HitWallCheck(this,surface_zIntercept,surface_dz_dx,tt)
            
            if this.x > ( this.z - surface_zIntercept ) / surface_dz_dx;
                
                this.amu = 1e200;
                this.vx = 0;
                this.vy = 0;
                this.vz = 0;
                this.hitWall = 1;
                
            end
            
        end
        
        function OutOfDomainCheck(this,xMinV,xMaxV,yMinV,yMaxV,zMinV,zMaxV)
            
            if this.x > xMaxV || this.x < xMinV ...
                    || this.y > yMaxV || this.y < yMinV ...
                    || this.z > zMaxV || this.z < zMinV;
                
                this.amu = 1e200;
                this.vx = 0;
                this.vy = 0;
                this.vz = 0;
                this.leftVolume = 1;
                
            end
        end
        
        function UpdatePrevious(this)
            this.xPrevious = this.x;
            this.yPrevious = this.y;
            this.zPrevious = this.z;
            
            this.vxPrevious = this.vx;
            this.vyPrevious = this.vy;
            this.vzPrevious = this.vz;
        end
        
        function PerpDistanceToSurface(this,surface_dz_dx,surface_zIntercept)
            if this.hitWall == 0 && this.leftVolume ==0
            this.perpDistanceToSurface = ( -surface_dz_dx*this.x + this.z + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
            end
        end
        
        function Energy(this)
                        MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            this.energy = 0.5*184*MI*(this.vxPrevious^2 +this.vyPrevious^2 + this.vzPrevious^2)/Q;
        end
        
       
       
    end
end

