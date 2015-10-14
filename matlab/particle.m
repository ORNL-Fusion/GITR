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
        function status = move (p, t, dt, E, B, xyz)
            
            
            status = 0;
            
            y0 = [p.x p.y p.z p.vx p.vy p.vz];
            tSpan = [0,t];
            
            options = odeset('InitialStep',dt,'MaxStep',dt);
                        
            [t_,y] = ode45(@(t_,y0) myode(t_,y0,p,B,xyz),tSpan,y0,options);
            
        end
    end 

end

