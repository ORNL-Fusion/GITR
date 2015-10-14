classdef particle
    
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
        function status = move (p, dt)
            
            constants
            
            status = 0;
            
            r0 = [p.x p.y p.z]';
            v0 = [p.vx p.vy p.vz]';
            y0 = [r0; v0];
            tSpan = [0,dt];
            this_q = p.Z * q;
            E = [0 0 0]';
            B = [0 0 0]';
            m = p.amu * mi;
            
            f = @(t,y) [y(4:6); (this_q/m) * (E + cross(y(4:6),B) )];
            
            [t,y] = ode45(f,tSpan,y0);
            
            p.x = y(1);
            p.y = y(2);
            p.z = y(3);
            p.vx = y(4);
            p.vy = y(5);
            p.vz = y(6);
            
        end
    end 

end

