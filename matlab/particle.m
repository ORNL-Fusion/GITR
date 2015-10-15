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
        function [T Y] = move (p, end_t, dt, E, B, xyz)
            
            
            status = 0;
            
            IC = [p.x p.y p.z p.vx p.vy p.vz]';
            tSpan = [0,end_t];
            
            options = odeset('InitialStep',dt,'MaxStep',dt);
                        
            [T Y] = ode45(@(t,y) myode(t,y,p,E,B,xyz),tSpan,IC,options);
      
        end
    end 

end

