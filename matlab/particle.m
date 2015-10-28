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



        function boris(part,xV,yV,zV,Bx,By,Bz,BMag,Ex,Ey,Ez,dt)
            q = 1.602e-19;
            amu_mass = 1.66e-27;
            BxPart = interpn(xV,yV,zV,Bx,part.x,part.y,part.z);
            ByPart = interpn(xV,yV,zV,By,part.x,part.y,part.z);
            BzPart = interpn(xV,yV,zV,Bz,part.x,part.y,part.z);
            BMagPart = sqrt( BxPart.^2 + ByPart.^2 + BzPart.^2 );
            
            ExPart = interpn(xV,yV,zV,Ex,part.x,part.y,part.z);
            EyPart = interpn(xV,yV,zV,Ey,part.x,part.y,part.z);
            EzPart = interpn(xV,yV,zV,Ez,part.x,part.y,part.z);
            
            %                 %Constants used in Boris method Lorentz Integrator
         q_prime = part.Z*q/(part.amu*amu_mass)*dt/2;
        coeff = 2*q_prime/(1+(q_prime*BMagPart).^2);

        %Boris Method Lorentz Integrator
        v = [part.vx part.vy part.vz];
        E = [ExPart EyPart EzPart];
        B = [BxPart ByPart BzPart];
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
            q_m = q/m;
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
    
    end
end

