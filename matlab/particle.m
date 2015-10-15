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
    end
    
    methods
        function move(part,xV,yV,zV,Bx,By,Bz,BMag,Ex,Ey,Ez,dt)

            BxPart = interpn(xV,yV,zV,Bx,part.x,part.y,part.z)
            ByPart = interpn(xV,yV,zV,By,part.x,part.y,part.z)
            BzPart = interpn(xV,yV,zV,Bz,part.x,part.y,part.z)
         f = @(t,y) [[part.vx part.vy part.vz]'; part.Z/part.amu*cross([part.vx part.vy part.vz],[BxPart ByPart BzPart])'];
         [t,y] = ode45(f,[0 dt], [part.x; part.y; part.z; part.vx; part.vy; part.vz]);
t
size(t)
disp('y')
y(:,1:3)
size(y)

part.x = y(end,1)
part.y = y(end,2)
part.z = y(end,3)
part.vx = y(end,4);
part.vy = y(end,5);
part.vz = y(end,6)

plot3(y(:,1),y(:,2),y(:,3))
stop
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
         q_prime = part.Z*q/(part.amu*amu_mass)*dt/2
        coeff = 2*q_prime/(1+(q_prime*BMagPart).^2)

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
    

        
    end
    
end

