classdef particle < handle
    
    properties
        x
        y
        z
        vx
        vy
        vz
        Z
        amu
        hitWall
    end
    
    methods

        function [T Y] = move (p, end_t, dt, E, B, xyz)
            
            
            status = 0;
            
            IC = [p.x p.y p.z p.vx p.vy p.vz]';
            tSpan = [0,end_t];
            
            options = odeset('InitialStep',dt,'MaxStep',dt);
                        
            [T Y] = ode45(@(t,y) myode(t,y,p,E,B,xyz),tSpan,IC);%,options);
            
            p.x = Y(end,1); % this doesn't seem to work, i.e., cannot modify p outside this scope :(
            p.y = Y(end,2);
            p.z = Y(end,3);
            
            p.vx = Y(end,4);
            p.vy = Y(end,5);
            p.vz = Y(end,6);
      
        end



        function [E_last B_last] = boris(part,xV,yV,zV,Bx,By,Bz,BMag,Ex,Ey,Ez,dt,E_last,B_last)
            q = 1.602e-19;
            amu_mass = 1.66e-27;
            B(1) = interpn(xV,yV,zV,Bx,part.x,part.y,part.z);
            B(2) = interpn(xV,yV,zV,By,part.x,part.y,part.z);
            B(3) = interpn(xV,yV,zV,Bz,part.x,part.y,part.z);
            
            
            E(1) = interpn(xV,yV,zV,Ex,part.x,part.y,part.z);
            E(2) = interpn(xV,yV,zV,Ey,part.x,part.y,part.z);
            E(3) = interpn(xV,yV,zV,Ez,part.x,part.y,part.z);

            
            
            if any([isnan(E) isnan(B)]) == 1
                E = E_last;
                B = B_last;
            else
                E_last = E;
                B_last = B;
            end
            
            BMagPart =norm(B);
            %                 %Constants used in Boris method Lorentz Integrator
            q_prime = part.Z*q/(part.amu*amu_mass)*dt/2;
            coeff = 2*q_prime/(1+(q_prime*BMagPart).^2);
            
            %Boris Method Lorentz Integrator
            v = [part.vx part.vy part.vz];
            r = [part.x part.y part.z];
            
            v_minus = v + q_prime*E;
            
            v = v_minus + q_prime*cross(v_minus,B);
            
            v = v_minus + coeff*cross(v,B);
            
            v = v + q_prime*E;
            
            r = r + v*dt;

            part.x =r(1);
            part.y =r(2);
            part.z =r(3);
            part.vx = v(1);
            part.vy = v(2);
            part.vz = v(3);
        end
        
        function [T Yr] = rk4(part,xV,yV,zV,Bx,By,Bz,BMag,Ex,Ey,Ez,end_t,dt,steps)
            q = 1.602e-19;
            amu_mass = 1.66e-27;
            m = amu_mass*part.amu;
            q_m = part.Z*q/m;
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
                
            function [E B] = field_interp(r,xV,yV,zV,Bx,By,Bz,Ex,Ey,Ez)
                B(1) = interpn(xV,yV,zV,Bx,r(1),r(2),r(3));
                B(2) = interpn(xV,yV,zV,By,r(1),r(2),r(3));
                B(3) = interpn(xV,yV,zV,Bz,r(1),r(2),r(3));
                
                
                E(1) = interpn(xV,yV,zV,Ex,r(1),r(2),r(3));
                E(2) = interpn(xV,yV,zV,Ey,r(1),r(2),r(3));
                E(3) = interpn(xV,yV,zV,Ez,r(1),r(2),r(3));
            end
        
        end
        function  ionization(p,Bfield,Efield,xyz,dt,temp_eV,density,RateCoeff,Te,s3)
            
            amu_mass = 1.66e-27;
            q = 1.602e-19;
            %T = 0.5*p.amu*amu_mass*(p.vx^2 +p.vy^2+ p.vz^2)/q;
            
            minT = min(Te(:))
            T=interpn(xyz.x,xyz.y,xyz.z,temp_eV(:,:,:,1),p.x,p.y,p.z)
           if T > minT
            n=interpn(xyz.x,xyz.y,xyz.z,density(:,:,:,1),p.x,p.y,p.z)
            
            
            Coeff = interp1(Te(:,p.Z),RateCoeff(:,p.Z),T)
            
            if ( isnan(Coeff) )
                error('Ionization interpolation out of range')
            end
            
            tion = 1/(Coeff*n)
            P1 = 1-exp(-dt/tion)
            
            r1=rand(s3)
            if r1 <= P1
                p.Z = p.Z+1
            end
           end
        end
        
        function recombination(p,Bfield,Efield,xyz,dt,temp_eV,density,RecombRateCoeff,DTEVD,DDENSD,s4)



            amu_mass = 1.66e-27;
            q = 1.602e-19;
            %T = 0.5*p.amu*amu_mass*(p.vx^2 +p.vy^2+ p.vz^2)/q;
            
            minT = min(DTEVD)
            minN = min(DDENSD)
            T=interpn(xyz.x,xyz.y,xyz.z,temp_eV(:,:,:,1),p.x,p.y,p.z)
            n=interpn(xyz.x,xyz.y,xyz.z,density(:,:,:,1),p.x,p.y,p.z)
            T = log10(T)
            n = log10(n)
           if (T > minT) && (n > minN) && p.Z >0
            





Coeff = interpn(DTEVD,DDENSD,RecombRateCoeff(:,:,p.Z),T,n);

tion = 1/(Coeff*10^n)
P1 = exp(dt/tion)


 if rand(s4) >= P1
    p.Z = p.Z-1
 end   
           end
        end
        
        function slow(part,xyz,Bx,By,Bz,BMag,Ex,Ey,Ez,temp_eV,density,dt,nS,amu,Z,s5,s6,s7)
            e0 = 8.85e-12;
            mi = 1.66e-27;
            q = 1.602e-19;
            
            E = [0 0 0];
            B = [0 0 0];
            v = [part.vx part.vy part.vz];
            r = [part.x part.y part.z];
            
            T=interpn(xyz.x,xyz.y,xyz.z,temp_eV(:,:,:,1),part.x,part.y,part.z)
            n=interpn(xyz.x,xyz.y,xyz.z,density(:,:,:,1),part.x,part.y,part.z)
            
            flow_v = [0 0 0]; %This will need to be calculated or interpolated in the future
            v_relative = v - flow_v;
            v_norm = norm(v_relative);
            [E B] = field_interp(r,xyz,Bx,By,Bz,Ex,Ey,Ez);
            nu_totals = [0 0 0];
            
            z = part.Z;
            m = part.amu*mi;
            for j=1:nS
              
                zbackground = Z(j);
                
                mbackground = amu(j)*mi
                
                lam_d = sqrt(e0*T/(n*zbackground^2*q));%only one q in order to convert to J
                lam = 4*pi*n*lam_d^3;
                gam = q^4*z^2*zbackground^2*log(lam)/(m*m*4*pi*e0*e0);
                a = mbackground/(2*T*q); %q is just to convert units - no z needed
                
                x = v_norm^2*a
                psi_prime = 2*sqrt(x/pi)*exp(-x)
                psi_psiprime = erf(sqrt(x))
                psi = psi_psiprime - psi_prime
                nu_0 = gam*n/v_norm^3
                nu_s = -(1+m/mbackground)*psi*nu_0
                nu_pitchangle = 2*(psi_psiprime - psi/(2*x))*nu_0
                nu_par = psi/x*nu_0
                nu_E = 2*(m/mbackground*psi - psi_prime)*nu_0
                nu_E =2*(m/mbackground*psi)*nu_0
                nu_totals(1) = nu_totals(1)+nu_s;
                nu_totals(2) = nu_totals(2)+nu_pitchangle;
                nu_totals(3) = nu_totals(3) + nu_E;
            end
            nu_totals
            [e1 e2 e3] = direction(B, v_relative);
            dv_relax = [0 0 0];
            dv_relax = dv_relax + e1*(nu_totals(1)*v_norm*dt + sqrt(nu_totals(3)*v_norm^2*dt));
            %nu_s*norm(v_relative)*dt + sqrt(nu_E*norm(v_relative)^2*dt)
            %norm(dv_relax)
            dv_relax = dv_relax + e2*sqrt(nu_totals(2)/2*v_norm^2*dt);
            dv_relax = dv_relax + e3*sqrt(nu_totals(2)/2*v_norm^2*dt)
            function [E B] = field_interp(r,xyz,Bx,By,Bz,Ex,Ey,Ez)
                B(1) = interpn(xyz.x,xyz.y,xyz.z,Bx,r(1),r(2),r(3));
                B(2) = interpn(xyz.x,xyz.y,xyz.z,By,r(1),r(2),r(3));
                B(3) = interpn(xyz.x,xyz.y,xyz.z,Bz,r(1),r(2),r(3));
                
                
                E(1) = interpn(xyz.x,xyz.y,xyz.z,Ex,r(1),r(2),r(3));
                E(2) = interpn(xyz.x,xyz.y,xyz.z,Ey,r(1),r(2),r(3));
                E(3) = interpn(xyz.x,xyz.y,xyz.z,Ez,r(1),r(2),r(3));
            end
            
            function [e1 e2 e3] = direction(B, v_relative)
                B_unit = B/norm(B);
                g = v_relative;
                e3 = g/norm(g);
                
                s1 = dot(e3,B_unit);
                s2 = sqrt(1-s1^2);
                
                e1 = 1/s2*(s1*e3 - B_unit);
                e2 = -1/s2*cross(e3,B_unit);
            end
            end
        
    end
end

